"""All mechanical and thermal properties assume that x- and y-axes are in the
   0- and 90-degree directions, respectively"""

import numpy as np
import math
import csv


class Lamina:
    """Individual lamina"""

    def __init__(self, plyID, thickness, theta, E11, E22, nu12, G12, a11, a22,
                 b11, b22):
        self.plyID = plyID
        self.thickness = thickness
        self.z = 0
        self.theta = theta  # theta must be in radians
        self.E11 = E11
        self.E22 = E22
        self.nu12 = nu12
        self.G12 = G12
        self.a11 = a11
        self.a22 = a22
        self.b11 = b11
        self.b22 = b22

        # create stiffness matrix and temp values
        self.matQ = np.zeros((3, 3))

        nu21 = nu12 * E22 / E11
        q11 = E11 / (1 - nu12 * nu21)
        q12 = nu12 * E22 / (1 - nu12 * nu21)
        q22 = E22 / (1 - nu12 * nu21)
        q66 = G12

        cost = math.cos(self.theta)
        sint = math.sin(self.theta)

        # build the transformed stiffness matrix for the ply
        self.matQ[0][0] = q11 * cost ** 4 + 2 * (q12 + 2 * q66) * cost ** 2 * \
                          sint ** 2 + q22 * sint ** 4
        self.matQ[0][1] = q12 * (cost ** 4 + sint ** 4) + (q11 + q22 - 4 *
                          q66) * cost ** 2 * sint ** 2
        self.matQ[0][2] = (q11 - q12 - 2 * q66) * cost ** 3 * sint - \
                          (q22 - q12 - 2 * q66) * cost * sint ** 3
        self.matQ[1][0] = self.matQ[0][1]
        self.matQ[1][1] = q11 * sint ** 4 + 2 * (q12 + 2 * q66) * sint ** 2 * \
                          cost ** 2 + q22 * cost ** 4
        self.matQ[1][2] = (q11 - q12 - 2 * q66) * cost * sint ** 3 - \
                          (q22 - q12 - 2 * q66) * cost ** 3 * sint
        self.matQ[2][0] = self.matQ[0][2]
        self.matQ[2][1] = self.matQ[1][2]
        self.matQ[2][2] = (q11 + q22 - 2 * q12 - 2 * q66) * cost ** 2 * \
                          sint ** 2 + q66 * (sint ** 4 + cost ** 4)

        self._calcABD()

    def _calcABD(self):
        """Calculate the A, B , and D matrices for the lamina within it's
        laminate"""

        # Calculate the upper and lower planes of the laminate
        zk = self.z + self.thickness / 2
        zk1 = self.z - self.thickness / 2

        self.matA = self.matQ * (zk - zk1)
        self.matB = (1 / 2) * self.matQ * (zk**2 - zk1**2)
        self.matD = (1 / 3) * self.matQ * (zk**3 - zk1**3)


class Laminate(list):
    """Laminate list object made up of multiple plies of lamina"""

    laminatesCount = 0

    def __init__(self):
        """Initialize with a list of Lamina objects"""

        # preserve list __init__() and append with additional laminate
        # properties
        super().__init__()

        Laminate.laminatesCount += 1
        self.laminateID = Laminate.laminatesCount

        # determine laminate properties for non-empty Laminate object
        if len(self) > 0:
            self.__lamUpdate()

    def __lamUpdate(self):
        """Calculates the ABD matrix and determines derivative effective
        laminate properties"""

        self.thickness = 0  # reset thickness
        self.layup = []
        self.thickness = 0
        self.matA = np.zeros((3, 3))
        self.matB = np.zeros((3, 3))
        self.matD = np.zeros((3, 3))

        # calculate z for each ply with respect to the bottom plane and
        # laminate thickness
        for ply in self:
            ply.z = self.thickness + ply.thickness / 2
            self.thickness += ply.thickness

        # recalculate z with respect to the laminate mid-plane
        for ply in self:
            ply.z -= self.thickness / 2
            ply._calcABD()

        # determine layup
        for ply in self:
            self.layup.append(math.degrees(ply.theta))

        # calculate A, B, and D matrices
        for ply in self:  # sequence ^
            self.matA += ply.matA
            self.matB += ply.matB
            self.matD += ply.matD

        # rebuild ABD matrix
        self.matABD = np.r_[np.c_[self.matA, self.matB],
                            np.c_[self.matB, self.matD]]

        # recalculate mechanical properties
        # axial properties
        invA = np.linagl.inv(self.matA)
        self.Ex = 1 / (invA[0][0] * self.thickness)
        self.Ey = 1 / (invA[1][1] * self.thickness)
        self.Gxy = 1 / (invA[2][2] * self.thickness)
        self.nuxy = - invA[0][1] / invA[1][1]
        self.nuyx = - invA[0][1] / invA[0][0]

        # flexural properties
        invD = np.linalg.inv(self.matD)
        self.Efx = 12 / (invD[0][0] * self.thickness**3)
        self.Efy = 12 / (invD[1][1] * self.thickness**3)
        self.Gfxy = 12 / (invD[2][2] * self.thickness**3)
        self.nufxy = - invD[0][1] / invD[1][1]
        self.nufyx = - invD[0][1] / invD[0][0]

    def append(self, newPly):
        """Extend list.append() to update laminate properties on addition of
           new ply. Note: added plies are assumed to be placed on TOP SURFACE
           of existing laminate."""

        # preserve list.append() and extend to update effective laminate
        # properties when a new ply is added
        super().append(newPly)
        self.__lamUpdate()

    def __repr__(self):
        return "Laminate ID: %s\n" \
               "Number of plies: %s\n" \
               "Thickness: %s\n" \
               "Layup: %s\n" \
               "Extensional Properties\n" \
               "    Ex: %s\n" \
               "    Ey: %s\n" \
               "    Gxy: %s\n" \
               "    nuxy: %s\n" \
               "    nuyx: %s\n" \
               "Flexural Properties\n" \
               "    Efx: %s\n" \
               "    Efy: %s\n" \
               "    Gfxy: %s\n" \
               "    nufxy: %s\n" \
               "    nufyx: %s\n" \
               % (self.laminateID,
                  len(self),
                  self.thickness,
                  self.layup,
                  self.Ex,
                  self.Ey,
                  self.Gxy,
                  self.nuxy,
                  self.nuyx,
                  self.Efx,
                  self.Efy,
                  self.Gfxy,
                  self.nufxy,
                  self.nufyx,)


def LamPropFromCSV(inputFile):
    """Determines laminate properties from input CSV file."""

    with open(inputFile, 'r') as csvLamina:
        rawLines = csv.DictReader(csvLamina)  # create dict object from csv

        # sort by plyID - assumed order is from the bottom upward
        sortedLines = sorted(rawLines,
                             key=lambda sortField: int(sortField['plyID']))

        L = Laminate()
        for ply in sortedLines:
            L.append(Lamina(int(ply['plyID']),
                            float(ply['thick']),
                            math.radians(float(ply['theta'])),
                            float(ply['E11']),
                            float(ply['E22']),
                            float(ply['nu12']),
                            float(ply['G12']),
                            float(ply['a11']),
                            float(ply['a22']),
                            float(ply['b11']),
                            float(ply['b22'])))
