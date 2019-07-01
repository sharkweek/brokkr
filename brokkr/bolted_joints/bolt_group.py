"""Take user-supplied fastener information fastener information (i.e. location,
diameter, etc.) and applied loads for a 3D fastener group and solve for the
loads applied to each fastener.

Version:
    Python 3.x

Required modules:
    - np
    - math
    - csv

Inputs:
    On running the script, the user is requested to input the file path for
    input files conatining the the necessary fastener and loading information.

    Fastener and applied loads information must be in a .CSV file. The
    fastener information .CSV must contain the following fields and headers:
        ID (int): fastener ID (must be unique integers)
        x, y, z (float): X, Y, and Z coordinates for each fastener
        dia (float): fastener diameter
        E (float): Young's Modulus
        G (float): shear modulus
        l (float): grip length

    The applied loads .CSV must contain the following fields and headers
    (in order):
        ID (int) - load ID (each ID must be unique)
        x, y, z (float): X, Y, and Z coordinates for each applied load
        fx, fy, fz (float): X, Y, and Z force components for each applied load
        mx, my, mz (float): X, Y, and Z moment components for each applied
            load

Assumptions:
    - rigid body mechanics (see TODO)
    - unit system remains consistent (i.e. US or SI)

Returns:
    fastener_loads.csv (file): CSV containing all the force vectors for each
        fastener

TODO:
    - Script currently accepts fastener moduli, but conservatively assumes
      rigid body mechanics. Functionality needs to be added to calculate loads
      based on fastener stiffness
    - Add 'Totals' row to end of output file
    - report values smaller than 1e-9 as zeros

"""

import numpy as np
from numpy import array, zeros
import math
import csv
from .fastener import Fastener
from .load import Load
from .loadset import LoadSet


class FastenerGroup(dict):
    """A fastener group."""

    # TODO: follow Niu Stress analysis chapter 9.4

    def __init__(self, *fasteners, load=None):
        """Initialize the instance."""

        # make sure each fastener is a Fastener class
        for each in fasteners:
            self.__is_fastener(each)

        super().__init__({i + 1: j for i, j in enumuerate(fasteners)})

        self.load = load

        self.__update()

    def __delitem__(self, key):
        """Delete a fastener at the specified index."""

        super().__delitem__(key)
        self.__update()

    def __is_fastener(f):
        """Check if object is a Fastener."""

        if type(f) != Fastener:
            raise TypeError("`FastenerGroup`s may contain only Fasteners")
        else:
            return True

    def __setitem__(self, key, new_fastener):
        """Set an item at a specified index to new_fastener."""

        self.__is_fastener(new_fastener)
        self[key] = new_fastener
        self.__update()

    def clear(self):
        """Clear the FastenerGroup of all Fasteners."""

        super().clear()
        self.__update()

    @property
    def cg(self):
        """The center of gravity of the fastener group.

        The x-component of the CG is calculated

            X_bar = sum(weight * x_i) / sum(weight)

        The weight factor 'weight' is typically a function of stiffness (i.e.
        area multiplied by modulus), but in this case it is a function of
        Fastener.wt.
        """

        cg = zeros((3, 1))
        sum_product = zeros((3, 1))
        sum_weight = 0

        for fstnr in self:
            # Calculate the sum of the products
            for i, component in enumerate(self[fstnr].xyz):
                sum_product[i] += component * self[fstnr].wt

            # Calculate the sum of the areas
            sum_weight += self[fstnr].wt

        # Divide sum of products by sum of areas
        for i, product in enumerate(sum_product):
            cg[i] = product / sum_weight

        return cg

    @property
    def moi(self):
        """Return the second moment of inertia for the group about each axis.

        Each value in the returned vector corresponds to the app.
        """

        pass

    def __update(self):
        """Update bolt group to distribute loads."""

        # Make sure loads have been assigned to group
        if type(self.load) == Load:
            self.load = LoadSet(self.load)
        elif type(self.load) != LoadSet:
            raise TypeError("Applied load must be a Load or LoadSet")

        # Begin Calculations
        _cg = self.cg  # calculate the cg once to save computation time
        _appLoad = self.load.totalForce
        _appMoment = self.load.totalMoment

        coef_mat = np.zeros((len(self) * 3, len(self) * 3))  # coeff matrix
        soln_mat = np.zeros(len(self) * 3)  # solution matrix

        cSet = [[i, i+1, i+2] for i in range(0, 3 * len(self), 3)]
        rSet = [[i+6, i+7, i+8] for i in range(0, 3 * (len(self) - 2), 3)]

        for i, j in enumerate(cSet):
            # i = column fastener ID
            # j = column fastener set
            # Mx = yFz - zFy
            # My = zFx - xFz
            # Mz = xFy - yFx

            Fx = j[0]
            Fy = j[1]
            Fz = j[2]

            # fill in first three rows
            coef_mat[0][Fx] = 1  # sum of Fx
            coef_mat[1][Fy] = 1  # sum of Fy
            coef_mat[2][Fz] = 1  # sum of Fz

            # fill in fourth row (sum of Mx at CG)
            coef_mat[3][Fy] = -(F[i].xyz[2] - _cg[2])  # -zFy
            coef_mat[3][Fz] = +(F[i].xyz[1] - _cg[1])  # +yFz

            # fill in fifth row (sum of My at CG)
            coef_mat[4][Fx] = +(F[i].xyz[2] - _cg[2])  # +zFx
            coef_mat[4][Fz] = -(F[i].xyz[0] - _cg[0])  # -xFz

            # fill in sixth row (sum of Mz at CG)
            coef_mat[5][Fx] = -(F[i].xyz[1] - _cg[1])  # -yFx
            coef_mat[5][Fy] = +(F[i].xyz[0] - _cg[0])  # +xFy

            for u, w in enumerate(rSet):
                # u = row fastener ID
                # w = row fastener set

                rX = w[0]
                rY = w[1]
                rZ = w[2]

                coef_mat[rX][Fy] = -(F[i].xyz[2] - F[u].xyz[2])  # -zFy
                coef_mat[rX][Fz] = +(F[i].xyz[1] - F[u].xyz[1])  # +yFz

                coef_mat[rY][Fx] = +(F[i].xyz[2] - F[u].xyz[2])  # +zFx
                coef_mat[rY][Fz] = -(F[i].xyz[0] - F[u].xyz[0])  # -xFz

                coef_mat[rZ][Fx] = -(F[i].xyz[1] - F[u].xyz[1])  # -yFx
                coef_mat[rZ][Fy] = +(F[i].xyz[0] - F[u].xyz[0])  # +xFy

        # fill in the solution matrix (soln_mat)
        for i in range(3):
            soln_mat[i] = -_netLoad.force[i]
            soln_mat[i+3] = -_netLoad.moment[i]

        # fill in the remaining rows
        for i, j in enumerate(rSet):
            # i = fastener
            # j = row

            rX = j[0]
            rY = j[1]
            rZ = j[2]

            # Mx = (y_cg - y_i)F_znet - (z_cg - z_i)F_ynet + M_xnet
            soln_mat[rX] = - ((_cg[1] - F[i].xyz[1]) * _netLoad.force[2]
                          - (_cg[2] - F[i].xyz[2]) * _netLoad.force[1]
                          + _netLoad.moment[0])

            # My = (z_cg - z_i)F_xnet - (x_cg - x_i)F_znet + M_ynet
            soln_mat[rY] = -((_cg[2] - F[i].xyz[2]) * _netLoad.force[0]
                           - (_cg[0] - F[i].xyz[0]) * _netLoad.force[2]
                           + _netLoad.moment[1])

            # Mz = (x_cg - x_i)F_ynet - (y_cg - y_i)F_xnet + M_znet
            soln_mat[rZ] = -((_cg[0] - F[i].xyz[0]) * _netLoad.force[1]
                         - (_cg[1] - F[i].xyz[1]) * _netLoad.force[0]
                         + _netLoad.moment[2])

        # Solve system of equations
        matSol = np.linalg.lstsq(coef_mat, soln_mat)[0]

        # Add resulting fastener loads to fastener objects
        for i, j in enumerate(cSet):
            rX = j[0]
            rY = j[1]
            rZ = j[2]

            F[i].force[0] = matSol[rX]
            F[i].force[1] = matSol[rY]
            F[i].force[2] = matSol[rZ]


    def getfromCSV(self, fastPath):
        """Import list of fasteners from CSV"""

        with open(fastPath, "r") as inputFasteners:
            fastenerLines = csv.DictReader(inputFasteners)

            F = []
            # Create a list of Fastener objects from each line in the mapping file
            for fstnr in fastenerLines:
                # Parse each line for location, load, and diameter information
                fastID = int(fstnr['ID'])
                location = [float(fstnr['x']), float(fstnr['y']), float(fstnr['z'])]
                diameter = float(fstnr['dia'])
                E = float(fstnr['E'])
                G = float(fstnr['G'])
                l = float(fstnr['l'])

                # Replace each list item in F with Fastener objects
                F.append(Fastener(fastID, location, diameter, E, G, l))

    def writetoCSV(self, fileName):
        """Output resulting fastener loads to a CSV"""

        with open(fileName, 'w') as writeFile:
            writeFile.write("ID,Fx,Fy,Fz\n")
            for fstnr in F:
                writeFile.write(str(fstnr.ID))
                for i in fstnr.force:
                    writeFile.write(',' + str(i))
                writeFile.write('\n')
