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

import math
import numpy as np
import pandas as pd
import csv


class Fastener:
    """Fastener class containing an ID, location in space, and load vector."""

    def __init__(self,
                 ID,
                 xyz,
                 diameter,
                 E,
                 G,
                 length,
                 weighting=1,
                 axis=[0, 0, 0]):
        """Initialize the instance."""

        self.ID = ID
        self.xyz = xyz
        self.diameter = diameter
        self.E = E
        self.G = G
        self.length = length
        self.force = [0, 0, 0]
        self.weighting = weighting
        self.axis = axis

    def __repr__(self):
        """Return the "official" Fastener string representation."""

        return "Fastener ID: %s\nLocation: %s\nForce: %s\nDiameter: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.diameter)

    @property
    def area(self):
        """Cross-sectional area of the fastener."""

        return (self.diameter ** 2) * math.pi / 4

    @area.setter
    def area(self, a):
        """Determine the diameter if the area is set manually."""

        self.diameter = math.sqrt(4 * a / math.pi)


class Load:
    """Load class containing location, force, and moment vectors."""

    def __init__(self,
                 xyz=[0, 0, 0],
                 force=[0, 0, 0],
                 moment=[0, 0, 0]):
        """Initialize a Load instance."""

        self.xyz = xyz
        self.force = force
        self.moment = moment

    def __repr__(self):
        """Return the "official" Load string representation."""

        return "Load ID: %s\nLocation: %s\nForce: %s\nMoment: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.moment)


class LoadSet:
    """A set of loads."""

    def __init__(self, *loads):
        """Initilize a LoadSet instance."""

        # check that loads are Load objects
        for each in loads:
            self.__isLoad(each)
        self.__loads = list(loads)

        self.__sumLocation = [0, 0, 0]  # location about which to sum loads

        if len(self) > 0:
            self.__update()

    def __delitem__(self, key):
        """Delete a load with the specified index within the load set."""

        del self.__loads[key]
        self.__update()

    def __getitem__(self, key):
        """Return a load at the specified index."""

        return self.__loads[key]

    def __isLoad(load):
        """Check if object is a Load type."""

        if type(load) != Load:
            raise TypeError("LoadSets must be composed of Load objects")
        else:
            return True

    def __iter__(self):
        """Return the loads as iterator."""

        self.__counter = 0
        return iter(self.__loads)

    def __len__(self):
        """Return the number of loads in LoadSet."""

        return len(self.__loads)

    def __next__(self):
        """Iterate to the next load."""

        self.__counter += 1
        if self.__counter < len(self):
            return self.__loads[self.__counter]
        else:
            raise StopIteration

    def __setitem__(self, key, newLoad):
        """Overwrite a load at the specified index."""

        self.__isLoad(newLoad)  # Checks if new load is a Load type
        self.__loads[key] = newLoad
        self.__update()

    def __update(self):
        """Update totals based on Loads in set and __sumLocation."""

        # sum the forces
        self.__totalForce = [0, 0, 0]
        for load in self:
            self.__totalForce += load.force

        # sum of the moments
        self.__totalMoment = [0, 0, 0]
        for load in self:
            self.__totalMoment += load.moment

        # translate and add moments induced by forces
        # Mx = Sum(yFz - (zFy-CGy))
        # My = Sum(zFx - (xFz-CGz))
        # Mz = Sum(xFy - (yFx-CGx))
        for load in self:
            self.__totalMoment[0] += (load.force[1] * (self.sumLocation[2] -
                                                       load.xyz[2]) -
                                      load.force[2] * (self.sumLocation[1] -
                                                       load.xyz[1]))
            self.__totalMoment[1] += (load.force[2] * (self.sumLocation[0] -
                                                       load.xyz[0]) -
                                      load.force[0] * (self.sumLocation[2] -
                                                       load.xyz[2]))
            self.__totalMoment[2] += (load.force[0] * (self.sumLocation[1] -
                                                       load.xyz[1]) -
                                      load.force[1] * (self.sumLocation[0] -
                                                       load.xyz[0]))

    @property
    def totalForce(self):
        """Return __totalForce."""

        return self.__totalForce

    @totalForce.setter
    def totalForce(self, *args):
        """Raise error when totalForce is manually set."""

        raise AttributeError("""totalForce cannot be assigned manually.""")

    @property
    def totalMoment(self):
        """Return __totalMoment."""

        return self.__totalMoment

    @totalMoment.setter
    def totalMoment(self, *args):
        """Raise error when totalMoment is manually set."""

        raise AttributeError("""totalMoment cannot be assigned manually.""")

    def append(self, newLoad):
        """Add a new load to the load set."""

        self.__isLoad(newLoad)  # make sure newLoad is a Load object
        self.__loads.append(newLoad)
        self.__update()

    def get_from_csv(self, fileName, clear=False):
        """Import applied loads from a CSV."""

        if clear:
            self.clear()

        appliedLoad = pd.DataFrame()
        appliedLoad.read_csv(fileName, dtype=float)

        # check that all necessary columns exist
        if 'x' not in appliedLoad:
            raise KeyError("No column labeled 'x' in CSV.")
        elif 'y' not in appliedLoad:
            raise KeyError("No column labeled 'y' in CSV.")
        elif 'z' not in appliedLoad:
            raise KeyError("No column labeled 'z' in CSV.")
        elif 'fx' not in appliedLoad:
            raise KeyError("No column labeled 'fx' in CSV.")
        elif 'fy' not in appliedLoad:
            raise KeyError("No column labeled 'fy' in CSV.")
        elif 'fz' not in appliedLoad:
            raise KeyError("No column labeled 'fz' in CSV.")
        elif 'mx' not in appliedLoad:
            raise KeyError("No column labeled 'mx' in CSV.")
        elif 'my' not in appliedLoad:
            raise KeyError("No column labeled 'my' in CSV.")
        elif 'mz' not in appliedLoad:
            raise KeyError("No column labeled 'mz' in CSV.")

        for index, row in appliedLoad.iterrows():
            loc = [row['x'], row['y'], row['z']]
            f = [row['fx'], row['fy'], row['fz']]
            m = [row['mx'], row['my'], row['mz']]

            self.append(Load(loc, f, m))

    def insert(self, key, newLoad):
        """Insert new load at specified index."""

        self.__isLoad(newLoad)  # make sure newLoad is a Load object
        self.__loads.insert(key, newLoad)
        self.__update()

    def clear(self):
        """Clear the load set of loads."""

        self.__loads.clear()
        self.__update

    @property
    def sumLocation(self):
        """Return __sumLocation attribute.

        __sumLocation is kept private to ensure that, whenever it is changed,
        the set is updated and the total moment is accurate.
        """

        return self.__sumLocation

    @sumLocation.getter
    def sumLocation(self, point):
        """Update the total moment about a new location."""

        # make sure that the supplied point is a list of three numbers
        if type(point) != list and len(point) != 3:
            raise TypeError("Point must be a three-item list of numbers")
        else:
            for coord in point:
                if math.nan(coord):
                    raise TypeError("""Point must be a three-item list of
                                    numbers""")

        self.__sumLocation = point
        self.__update()


class FastenerGroup:
    """A fastener group."""

    def __init__(self, *fasteners):
        """Initialize the instance."""

        # make sure each fastener is a Fastener class
        if len(self) > 0:
            for each in fasteners:
                self.__isFastener(each)
            self.__fasteners = list(fasteners)

        self.__CG = [0, 0, 0]
        self.appliedLoad = None

        self.__update()

    def __delitem__(self, key):
        """Delete a fastener at the specified index."""

        del self.__fasteners[key]
        self.__update()

    def __len__(self):
        """Return number of fasteners in FastenerGroup."""

        return len(self.__fasteners)

    def __getitem__(self, key):
        """Return a Fastener at a specified index."""

        return self.__fasteners[key]

    def __update(self):
        """Update bolt group to distribute loads."""

        # Make sure loads have been specified
        if type(self.appliedLoad) == Load:
            self.appliedLoad = LoadSet(appliedLoad)
        elif type(self.appliedLoad) != LoadSet:
            raise TypeError("Applied load must be a Load or LoadSet")

        # Begin Calculations
        _cg = self.cg  # calculate the cg once to save computation time
        netLoad = self.appliedLoad.totalForce
        netMoment = self.appliedLoad.totalMoment

        matA = np.zeros((len(self) * 3, len(self) * 3))  # coeff matrix
        matB = np.zeros(len(self) * 3)  # solution matrix

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
            matA[0][Fx] = 1  # sum of Fx
            matA[1][Fy] = 1  # sum of Fy
            matA[2][Fz] = 1  # sum of Fz

            # fill in fourth row (sum of Mx at CG)
            matA[3][Fy] = -(F[i].xyz[2] - _cg[2])  # -zFy
            matA[3][Fz] = +(F[i].xyz[1] - _cg[1])  # +yFz

            # fill in fifth row (sum of My at CG)
            matA[4][Fx] = +(F[i].xyz[2] - _cg[2])  # +zFx
            matA[4][Fz] = -(F[i].xyz[0] - _cg[0])  # -xFz

            # fill in sixth row (sum of Mz at CG)
            matA[5][Fx] = -(F[i].xyz[1] - _cg[1])  # -yFx
            matA[5][Fy] = +(F[i].xyz[0] - _cg[0])  # +xFy

            for u, w in enumerate(rSet):
                # u = row fastener ID
                # w = row fastener set

                rX = w[0]
                rY = w[1]
                rZ = w[2]

                matA[rX][Fy] = -(F[i].xyz[2] - F[u].xyz[2])  # -zFy
                matA[rX][Fz] = +(F[i].xyz[1] - F[u].xyz[1])  # +yFz

                matA[rY][Fx] = +(F[i].xyz[2] - F[u].xyz[2])  # +zFx
                matA[rY][Fz] = -(F[i].xyz[0] - F[u].xyz[0])  # -xFz

                matA[rZ][Fx] = -(F[i].xyz[1] - F[u].xyz[1])  # -yFx
                matA[rZ][Fy] = +(F[i].xyz[0] - F[u].xyz[0])  # +xFy

        # fill in the solution matrix (matB)
        for i in range(3):
            matB[i] = -netLoad.force[i]
            matB[i+3] = -netLoad.moment[i]

        # fill in the remaining rows
        for i, j in enumerate(rSet):
            # i = fastener
            # j = row

            rX = j[0]
            rY = j[1]
            rZ = j[2]

            # Mx = (y_cg - y_i)F_znet - (z_cg - z_i)F_ynet + M_xnet
            matB[rX] = - ((_cg[1] - F[i].xyz[1]) * netLoad.force[2]
                          - (groupCG[2] - F[i].xyz[2]) * netLoad.force[1]
                          + netLoad.moment[0])

            # My = (z_cg - z_i)F_xnet - (x_cg - x_i)F_znet + M_ynet
            matB[rY] = -((_cg[2] - F[i].xyz[2]) * netLoad.force[0]
                         - (groupCG[0] - F[i].xyz[0]) * netLoad.force[2]
                         + netLoad.moment[1])

            # Mz = (x_cg - x_i)F_ynet - (y_cg - y_i)F_xnet + M_znet
            matB[rZ] = -((_cg[0] - F[i].xyz[0]) * netLoad.force[1]
                         - (groupCG[1] - F[i].xyz[1]) * netLoad.force[0]
                         + netLoad.moment[2])

        # Solve system of equations
        matSol = np.linalg.lstsq(matA, matB)[0]

        # Add resulting fastener loads to fastener objects
        for i, j in enumerate(cSet):
            rX = j[0]
            rY = j[1]
            rZ = j[2]

            F[i].force[0] = matSol[rX]
            F[i].force[1] = matSol[rY]
            F[i].force[2] = matSol[rZ]

    def __iter__(self):
        """Returns self as an iterator"""

        self.__counter = 0
        return iter(self.__fasteners)

    def __next__(self):
        """Iterates to the next item"""

        self.__counter += 1
        if self.__counter < len(self):
            return self.__fasteners[self.__counter]
        else:
            raise StopIteration

    def __setitem__(self, key, newFastener):
        """Sets an item at a specified index to newFastener"""

        self.__isFastener(newFastener)
        self.__fasteners[key] = newFastener
        self.__update()

    def append(self, newFastener):
        """Adds a newFastener to group"""

        self.__isFastener(newFastener)
        self.__fasteners.append(newFastener)
        self.__update()

    def clear(self):
        """Clears the FastenerGroup of all Fasteners"""

        self.__fasteners.clear()
        self.__update()

    @property
    def cg(self):
        """Finds the center of gravity of the fastener group"""
        _cg = [0, 0, 0]
        sumWeight = 0

        for fstnr in self:
            # Calculate the sum of the products
            for i, component in enumerate(fstnr.xyz):
                _cg[i] += component * fstnr.weighting

            # Calculate the sum of the areas
            sumWeight += fstnr.weighting

        # Divide sum of products by sum of areas
        for i, product in enumerate(_cg):
            _cg[i] = product / sumWeight

        return _cg

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

    def __isFastener(f):
        """Checks if object is a Fastener"""

        if type(f) != Fastener:
            raise TypeError("FastnerGroups may contain only Fasteners")
        else:
            return True

    def writetoCSV(self, fileName):
        """Output resulting fastener loads to a CSV"""

        with open(fileName, 'w') as writeFile:
            writeFile.write("ID,Fx,Fy,Fz\n")
            for fstnr in F:
                writeFile.write(str(fstnr.ID))
                for i in fstnr.force:
                    writeFile.write(',' + str(i))
                writeFile.write('\n')
