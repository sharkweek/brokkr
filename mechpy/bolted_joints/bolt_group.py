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
import csv


class Fastener:
    """Fastener class containing an ID, location in space, and load vector"""

    def __init__(self,
                 ID,
                 xyz,
                 diameter,
                 E,
                 G,
                 length,
                 axis=[0, 0, 0]):
        self.ID = ID
        self.xyz = xyz
        self.diameter = diameter
        self.E = E
        self.G = G
        self.length = length
        self.force = [0, 0, 0]
        self.axis = axis

    @property
    def area(self):
        """Cross-sectional area of the fastener"""

        return (self.diameter ** 2) * math.pi / 4

    @area.setter()
    def area(self, a):
        """Sets the diameter if the area is set by the user"""

        self.diameter = math.sqrt(4 * a / math.pi)

    def __repr__(self):
        return "Fastener ID: %s\nLocation: %s\nForce: %s\nDiameter: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.diameter)


class Load:
    """Load class containing location, force, and moment vectors"""

    def __init__(self,
                 xyz=[0, 0, 0],
                 force=[0, 0, 0],
                 moment=[0, 0, 0]):
        self.xyz = xyz
        self.force = force
        self.moment = moment

    def __repr__(self):
        return "Load ID: %s\nLocation: %s\nForce: %s\nMoment: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.moment)


class LoadSet:
    """A set of loads"""

    def isLoad(load):
        """Checks if object is a Load type"""

        if type(load) != Load:
            raise TypeError("LoadSets must be composed of Load objects")
        else:
            return True

    def __init__(self, loads):
        """Initilize the LoadSet"""

        # check that loads are given as a list
        if type(loads) != list:
            raise TypeError("Loads must be given in a list")
        else:
            for each in loads:
                self.isLoad(each)

        self.__loads = loads

        self.__update()

    def __len__(self):
        """Returns the number of loads in LoadSet"""

        return len(self.__loads)

    def __iter__(self):
        """Returns the loads as iterator"""

        self.__counter = 0
        return iter(self.__loads)

    def __next__(self):
        """Iterates to the next load"""

        self.__counter += 1
        if self.__counter < len(self):
            return self.__loads[self.__counter]
        else:
            raise StopIteration

    def clear(self):
        """Clears the load set of loads"""

        self.__loads.clear()
        self.__update

    def __getitem__(self, key):
        """Returns a load at the specified index"""

        return self.__loads[key]

    def append(self, newLoad):
        """Adds a new load to the load set"""

        self.__loads.append(newLoad)
        self.__update()

    def __delitem__(self, key):
        """Deletes a load with the specified index within the load set"""

        del self.__loads[key]
        self.__update()

    def __setitem__(self, key, newLoad):
        """Overwrites a load at the specified index"""

        self.isLoad(newLoad)  # Checks if new load is a Load type
        self.__loads[key] = newLoad
        self.__update()

    def __update(self):
        self.sumLoads():

    def sumLoads(self, point=[0, 0, 0]):
        """Determine the resultant load applied at specified point"""

        if type(point) != list and len(point) != 3:
            raise TypeError("Point must be a three-item list of numbers")
        else:
            for coord in point:
                if math.nan(coord):
                    raise TypeError("""Point must be a three-item list of
                                    numbers""")

        # sum the forces
        self.totalForce = [0, 0, 0]
        for load in self:
            self.totalForce += load.force

        # sum of the moments
        self.totalMoment = [0, 0, 0]
        for load in self:
            self.totalMoment += load.moment

        # translate and add moments induced by forces
        # Mx = Sum(yFz - (zFy-CGy))
        # My = Sum(zFx - (xFz-CGz))
        # Mz = Sum(xFy - (yFx-CGx))
        for load in self:
            self.totalMoment[0] += (load.force[1] * (point[2] - load.xyz[2]) -
                                    load.force[2] * (point[1] - load.xyz[1]))

            self.totalMoment[1] += (load.force[2] * (point[0] - load.xyz[0]) -
                                    load.force[0] * (point[2] - load.xyz[2]))

            self.totalMoment[2] += (load.force[0] * (point[1] - load.xyz[1]) -
                                    load.force[1] * (point[0] - load.xyz[0]))

    def fromCSV(filePath):
        """Import applied loads from a CSV"""

        with open(filePath, "r") as inputLoads:
            appliedload = csv.DictReader(inputLoads)

            self.clear()
            for ld in appliedload:
                loc = [float(ld['x']), float(ld['y']), float(ld['z'])]
                f = [float(ld['fx']), float(ld['fy']), float(ld['fz'])]
                m = [float(ld['mx']), float(ld['my']), float(ld['mz'])]

                self.append(Load(loc, f, m))


class FastenerGroup:
    """A fastener group"""

    def isFastener(f):
        """Checks if object is a Fastener"""

        if type(f) != Fastener:
            raise TypeError("FastnerGroups may contain only Fasteners")
        else:
            return True

    def __init__(self, fasteners):
        self.__fasteners = []

        if len(self) > 0:
            self.__groupUpdate_()

    def __len__(self):
        """Returns number of fasteners in FastenerGroup"""

        return len(self.__fasteners)

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

    def clear(self):
        """Clears the FastenerGroup of all Fasteners"""

        self.__fasteners.clear()
        self.__groupUpdate_()

    def append(self, newFastener):
        """Adds a newFastener to group"""

        self.isFastener(newFastener)
        self.__fasteners.append(newFastener)
        self.__groupUpdate_()

    def __getitem__(self, key):
        """Returns a Fastener at a specified index"""

        return self.__fasteners[key]

    def __delitem__(self, key):
        """Deletes a fastener at the specified index"""

        del self.__fasteners[key]
        self.__groupUpdate_()

    def __setitem__(self, key, newLoad)

    def fromCSV(self, fastPath):
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


    def findCG(self,fastenerList):
        """Finds the center of gravity of the fastener group"""
        CG = [0, 0, 0]
        sumArea = 0

        for fstnr in fastenerList:
            # Calculate the sum of the products
            for i, component in enumerate(fstnr.xyz):
                CG[i] += component * fstnr.area

            # Calculate the sum of the areas
            sumArea += fstnr.area

        # Divide sum of products by sum of areas
        for i, product in enumerate(CG):
            CG[i] = product / sumArea

        return CG

    def writetoCSV(self, fileName):
        """Output resulting fastener loads to a CSV"""

        with open(fileName, 'w') as writeFile:
            writeFile.write("ID,Fx,Fy,Fz\n")
            for fstnr in F:
                writeFile.write(str(fstnr.ID))
                for i in fstnr.force:
                    writeFile.write(',' + str(i))
                writeFile.write('\n')

    def __groupUpdate_(self):
        # Begin Calculations
        groupCG = findCG(F)
        netLoad = sumLoads(L, groupCG)
        matA = np.zeros((Fastener.count * 3, Fastener.count * 3))  # coeff matrix
        matB = np.zeros(Fastener.count * 3)  # solution matrix

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
        matSol = np.linalg.lstsq(matA, matB)[0]

        # Add resulting fastener loads to fastener objects
        for i, j in enumerate(cSet):
            rX = j[0]
            rY = j[1]
            rZ = j[2]

            F[i].force[0] = matSol[rX]
            F[i].force[1] = matSol[rY]
            F[i].force[2] = matSol[rZ]
