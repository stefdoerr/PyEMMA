__author__ = 'noe'


def ramachandran(X, order=2, ff='AMBER99'):
    """ Computes eigenvector basis functions as described in [5]_

    Parameters
    ---------
    X : ndarray(2,) or ndarray(T, 2)
        phi / psi angles for one residue at a given point in time, or for a full trajectory of length T
    order : int, default=2
        order of basis function famility. Here equal to the number of nontrivial MSM eigenvectors on phi/psi space
        to be used. The output will have order+1 columns
    ff : str, default='AMBER99'
        Force field to be used. Options: 'AMBER99' | 'CHARMM' | ...

    Returns
    -------
    Y : ndarray(T, order+1)
        The basis set values, evaluated for every trajectory time point and every function order.
        The first column is always the trivial 1-function.

    References
    ---------
    .. [5] Vitalini, F., Noe, F. and Keller, B. (2015):
        A basis set for peptides for the variational approach to conformational kinetics. (In review).

    """
    # TODO:
    # access table data corresponding to the chosen ff
    # read function values from the discrete phi/psi boxes corresponding to X (careful: must be vectorized over time - otherwise slow)
    # return Y
    pass

