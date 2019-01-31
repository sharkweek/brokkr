"""Load sets for bolted joint calculations."""

import math
from .load import Load


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

    @sumLocation.setter
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
