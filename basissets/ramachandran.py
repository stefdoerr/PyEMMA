from base import BasisSet

class RamachandranBasis(BasisSet):

    def __init__(self, restype, order=2, ff='AMBER99'):
        """Basis set on the phi/psi torsion space of a single amino acid.

            Parameters
            ---------
            restype : str
                Residue type, e.g. 'ALA'
            order : int, default=2
                order of basis function family. Here equal to the number of nontrivial MSM eigenvectors on phi/psi space
                to be used. The output will have order+1 columns
            ff : str, default='AMBER99'
                Force field to be used. Options: 'AMBER99' | 'CHARMM' | ...

            References
            ---------
            .. [5] Vitalini, F., Noe, F. and Keller, B. (2015):
                A basis set for peptides for the variational approach to conformational kinetics. (In review).

        """
        self.order = order
        self.ff = ff
        # TODO:
        # access table data corresponding to the chosen ff
        # read function values from the discrete phi/psi boxes corresponding to X (careful: must be vectorized over time - otherwise slow)


    def map(self, X):
        """
        Parameters
        ---------
        X : ndarray(T, 2)
            phi / psi angles for one residue in a  trajectory of length T

        Returns
        -------
        Y : ndarray(T, order+1)
            The basis set values, evaluated for every trajectory time point and every function order.
            The first column is always the trivial 1-function.

        """
        # evaluate phi/psi function (maybe use interpolation to get a smooth function)
        # return eigenvector values at phi/psi
        pass


class RamachandranProductBasis(BasisSet):

    def __init__(self, restypes, max_order=2, n_excite=-1, ff='AMBER99'):
        """Basis set on the phi/psi torsion space of a single amino acid.

        Parameters
        ---------
        restype : iterable of str
            Residue types along the chain, e.g. ['ALA', 'GLY', 'ALA', 'GLY', 'VAL']
        max_order : int, default=2
            maximum order of basis function family. Here equal to the number of nontrivial MSM eigenvectors
            on phi/psi space to be used. The output will have order+1 columns
        n_excite : int, default=-1
            number of Ramachandran basis functions allowed to be excited. The default value (-1) will allow
            all possible combinations (careful, this can be unmanageable). The minimal excitation is 1, where
            the excited states up to max_order are used for each amino acid in the chain.
        ff : str, default='AMBER99'
            Force field to be used. Options: 'AMBER99' | 'CHARMM' | ...

        References
        ---------
        .. [5] Vitalini, F., Noe, F. and Keller, B. (2015):
            A basis set for peptides for the variational approach to conformational kinetics. (In review).

        """
        self.restypes = restypes
        self.max_order = max_order
        self.n_excite = n_excite
        self.ff = ff
        # TODO:
        # build RamachandranBasis objects, one for each residues
        # build combinations according to max_order and n_excite


    def map(self, X):
        """
        Parameters
        ---------
        X : ndarray(T, 2R)
            phi / psi angles for each of R residues for a trajectory of length T

        Returns
        -------
        Y : ndarray(T, n)
            The basis set values, evaluated for every trajectory time point.
            Each of the n elements represents a product of the residue basis
            sets.

        """
        pass
