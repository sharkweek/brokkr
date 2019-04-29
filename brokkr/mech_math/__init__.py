"""Miscellaneous functions used in brokkr."""

import numpy as np
from .tensors import StrainTensor, StressTensor
from .vectors import ForceVector, MomentVector, DisplacementVector

__all__ = ['matrix_minor', 'ms', 'out_of_bounds']


def matrix_minor(matrix, indices):
    """Return the minor of a 2D ``numpy.ndarray``.

    Parameters
    ----------
    matrix : 2D numpy.ndarray
        the matrix from which to find a minor
    indices : tuple of ints
        a tuple of two values indicating the row and columns for the minor

    Returns
    -------
    numpy.ndarray
        the matrix minor

    Raises
    ------
    ValueError
        If the supplied ``indices`` is not a tuple of two integers or if
        ``matrix`` is not square

    """

    # error checking
    if len(indices) != 2:
        raise ValueError("`indices` must be a tuple of two integers")
    if matrix.ndim != 2 or len(matrix) != len(matrix[0]):
        raise ValueError("`matrix` must be a square 2D array")

    return np.delete(
        np.delete(matrix, indices[0], axis=0), indices[1], axis=1
    )


def ms(applied, allowed, knockdown=1):
    r"""Calculate the margin of safety.

    The margin of safety is calculated according to the following formula:

    .. math::
        MS = \frac{F_{\text{app}}}{k F_{alw}} - 1

    where

    ====================== ================
    :math:`MS`             margin of safety
    :math:`F_{\text{app}}` applied load
    :math:`F_{\text{alw}}` allowable
    ====================== ================

    .. note::
       The ``applied`` and ``allowed`` values must have the same units.

    Parameters
    ----------
    applied : float
        applied stress, load, or strain
    allowed : float
        the allowable strength
    knockdown : float
        allowable knockdown factor (``default=1``)

    Returns
    -------
    float
        the margin of safety

    Notes
    -----
    It is also possible to calculate a margin of safety using a failure index
    by substituting the failure index in for the ``applied`` value and setting
    the ``allowed`` value to ``1``

    """
    return (applied / (knockdown * allowed)) - 1


def out_of_bounds(val, mn, mx, condition):
    """Check if value satisfies boundary conditions.

    Parameters
    ----------
    val : float or int
        value to evaluate against boundary conditions
    mn, mx: float or int
        the minimum and maximum boundaries
    condition: {'g', 'ge', 'g-l', 'ge-l', 'g-le', 'ge-le', 'l', 'le'}
        the boundary condition to evaluate

    Returns
    -------
    bool
        True if ``val`` is outside boundaries, False if within boundaries

    Notes
    -----
    ``condition`` should may be any of the values defined for each of the
    boundary conditions described in the table below:

    =========== ==================
    Value       Boundary Condition
    =========== ==================
    ``'g'``     ``val > mn``
    ``'ge'``    ``val >= mn``
    ``'g-l'``   ``mn < val < mx``
    ``'ge-l'``  ``mn <= val < mx``
    ``'g-le'``  ``mn < val <= mx``
    ``'ge-le'`` ``mn <= val <= mx``
    ``'l'``     ``val < mx``
    ``'le'``    ``val <= mx``
    =========== ==================

    """

    return not {
        'g': lambda: val > mn and mx is None,
        'ge': lambda: val >= mn and mx is None,
        'g-l': lambda: mn < val < mx,
        'ge-l': lambda: mn <= val < mx,
        'g-le': lambda: mn < val <= mx,
        'ge-le': lambda: mn <= val <= mx,
        'l': lambda: mn is None and val < mx,
        'le': lambda: mn is None and val <= mx
    }.get(condition)()
