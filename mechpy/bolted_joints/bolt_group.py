"""Take user-supplied fastener information fastener information (i.e. location,
diameter, etc.) and applied loads for a 3D fastener group and solve for the
loads applied to each fastener.

Version:
    Python 3.x

Required modules:
    - numpy
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
    - Script currently accepts fastener moduli, but conservatively assumes rigid
      body mechanics. Functionality needs to be added to calculate loads based
      on fastener stiffness
    - Add 'Totals' row to end of output file
    - report values smaller than 1e-9 as zeros

"""

import math
import numpy
import csv


class Fastener(object):
    """Fastener class containing an ID, location in space, and load vector"""

    count = 0

    def __init__(self, ID, xyz, diameter, E, G, l):
        self.ID = ID
        self.xyz = xyz
        self.diameter = diameter
        self.E = E
        self.G = G
        self.l = l
        Fastener.count += 1
        self.area = (diameter ** 2) * math.pi / 4
        self.force = [0, 0, 0]

    def __repr__(self):
        return "Fastener ID: %s\nLocation: %s\nForce: %s\nDiameter: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.diameter)


class Load(object):
    """Load class containing location, force, and moment vectors"""

    count = 0

    def __init__(self, ID, xyz, force, moment):
        self.ID = ID
        self.xyz = xyz
        self.force = force
        self.moment = moment
        Load.count += 1

    def __repr__(self):
        return "Load ID: %s\nLocation: %s\nForce: %s\nMoment: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.moment)

    def loads_from_csv(self, loadsPath):

        # Request the mapping CSV file from the user. The format should be:
        # x,y,z,load_in_x,load_in_y,load_in_z,diameter
        fastPath = input("Fastener info *.csv file: ")
        loadsPath = input("Applied loads *.csv file: ")

        # Import applied Loads
        with open(loadsPath, "r") as inputLoads:
            appliedload = csv.DictReader(inputLoads)

            L = []
            for ld in appliedload:
                loadID = int(ld['ID'])
                loc = [float(ld['x']), float(ld['y']), float(ld['z'])]
                f = [float(ld['fx']), float(ld['fy']), float(ld['fz'])]
                m = [float(ld['mx']), float(ld['my']), float(ld['mz'])]

                L.append(Load(loadID, loc, f, m))

    def sumLoads(loadList, appLocation):
        """Determine the resultant load applied at the CG"""
        # sum the forces
        forces = [load.force for load in loadList]
        sumForce = [sum(f) for f in zip(*forces)]

        # sum of the applied moments
        moments = [load.moment for load in loadList]
        sumMoment = [sum(m) for m in zip(*moments)]

        # translate and add moments induced by forces
        # Mx = Sum(yFz - (zFy-CGy))
        # My = Sum(zFx - (xFz-CGz))
        # Mz = Sum(xFy - (yFx-CGx))
        for load in loadList:
            sumMoment[0] += load.force[1] * (appLocation[2] - load.xyz[2]) - \
                            load.force[2] * (appLocation[1] - load.xyz[1])

            sumMoment[1] += load.force[2] * (appLocation[0] - load.xyz[0]) - \
                            load.force[0] * (appLocation[2] - load.xyz[2])

            sumMoment[2] += load.force[0] * (appLocation[1] - load.xyz[1]) - \
                            load.force[1] * (appLocation[0] - load.xyz[0])

        sumLoadID = Load.count * 100 + 1

        return Load(sumLoadID, appLocation, sumForce, sumMoment)

class FastenerGroup(list):
    """A fastener group made up of multiple Fastener objects"""

    def __init__():
        """Extended the list.__init__()"""

        super().__init__()

        # determine laminate properties for non-empty Laminate object
        if len(self) > 0:
            self.__groupUpdate_()

    def append(self, newFastener):
        """Extended list.append() to update the fastener group properties on
        addition of new fastener."""

        # preserve list.append() and extend to update fastener group
        # properties when a new fastener is added
        super().append(newFastener)
        self.__groupUpdate_()

    def remove(self, delFastener):
        """Extended the list.remove() method to update the group when a
        fastener is removed"""

        super().remove(delFastener)
        self.__groupUpdate_()

    def insert(self, fastener):
        """Extended list.insert() to update laminate when a fastener is inserted"""

        super().insert(fastener)
        self.__groupUpdate_()

    def fasteners_from_csv(self, fastPath, overwrite=False):
        """Import to the fastener group from a CSV file"""

        # Clear the group from fasteners if desired
        if overwrite:
            self.clear()

        # Import fastener mapping file
        with open(fastPath, "r") as inputFasteners:
            fastenerLines = csv.DictReader(inputFasteners)

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

    @property
    def CG(self):
        """Finds the center of gravity of the fastener group"""

        CG = [0, 0, 0]
        sumArea = 0

        for fstnr in self:
            # Calculate the sum of the products
            for i, component in enumerate(fstnr.xyz):
                CG[i] += component * fstnr.area

            # Calculate the sum of the areas
            sumArea += fstnr.area

        # Divide sum of products by sum of areas
        for i, product in enumerate(CG):
            CG[i] = product / sumArea

        return CG

    def __groupUpdate_(self):
        """Updates the group properties"""

        # Begin Calculations
        groupCG = findCG(F)
        netLoad = sumLoads(L, groupCG)
        matA = numpy.zeros((Fastener.count * 3, Fastener.count * 3))  # coeff matrix
        matB = numpy.zeros(Fastener.count * 3)  # solution matrix

        cSet = [[i, i+1, i+2] for i in range(0, 3 * Fastener.count, 3)]
        rSet = [[i+6, i+7, i+8] for i in range(0, 3 * (Fastener.count - 2), 3)]

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
            matA[3][Fy] = -(F[i].xyz[2] - groupCG[2])  # -zFy
            matA[3][Fz] = +(F[i].xyz[1] - groupCG[1])  # +yFz

            # fill in fifth row (sum of My at CG)
            matA[4][Fx] = +(F[i].xyz[2] - groupCG[2])  # +zFx
            matA[4][Fz] = -(F[i].xyz[0] - groupCG[0])  # -xFz

            # fill in sixth row (sum of Mz at CG)
            matA[5][Fx] = -(F[i].xyz[1] - groupCG[1])  # -yFx
            matA[5][Fy] = +(F[i].xyz[0] - groupCG[0])  # +xFy

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
            matB[rX] = - ((groupCG[1] - F[i].xyz[1]) * netLoad.force[2]
                          - (groupCG[2] - F[i].xyz[2]) * netLoad.force[1]
                          + netLoad.moment[0])

            # My = (z_cg - z_i)F_xnet - (x_cg - x_i)F_znet + M_ynet
            matB[rY] = -((groupCG[2] - F[i].xyz[2]) * netLoad.force[0]
                       - (groupCG[0] - F[i].xyz[0]) * netLoad.force[2]
                       + netLoad.moment[1])

            # Mz = (x_cg - x_i)F_ynet - (y_cg - y_i)F_xnet + M_znet
            matB[rZ] = -((groupCG[0] - F[i].xyz[0]) * netLoad.force[1]
                       - (groupCG[1] - F[i].xyz[1]) * netLoad.force[0]
                       + netLoad.moment[2])

        # Solve system of equations
        matSol = numpy.linalg.lstsq(matA, matB)[0]

        # Add resulting fastener loads to fastener objects
        for i, j in enumerate(cSet):
            rX = j[0]
            rY = j[1]
            rZ = j[2]

            F[i].force[0] = matSol[rX]
            F[i].force[1] = matSol[rY]
            F[i].force[2] = matSol[rZ]

        # Write resulting fastener loads to output file
        with open('fastener_loads.csv', 'w') as writeFile:
            writeFile.write("ID,Fx,Fy,Fz\n")
            for fstnr in F:
                writeFile.write(str(fstnr.ID))
                for i in fstnr.force:
                    writeFile.write(',' + str(i))
                writeFile.write('\n')
