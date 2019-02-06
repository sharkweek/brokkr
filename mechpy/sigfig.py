"""A module for determining and rounding to significant figures."""

def sf(val):
    """
    Calculate the significant figures of a value

    Parameters
    ----------
    val : int or float
        value for which to calculate the significant figures
    
    """

    return len(str(val).strip