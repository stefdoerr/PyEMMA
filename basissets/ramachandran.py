from base import BasisSet
import numpy as np
import itertools as it

class RamachandranBasis(BasisSet):

    def __init__(self, restype, order=2, ff='AMBER99SB_ILDN'):
        """Basis set on the phi/psi torsion space of a single amino acid.

            Parameters
            ---------
            restype : str
                Residue type, single letter code, e.g. alanine = 'A'.
            order : int, default=2
                order of basis function family. Here equal to the number of nontrivial MSM eigenvectors on phi/psi space
                to be used. The output will have order+1 columns
            ff : str, default='AMBER99SB_ILDN'
                Force field to be used. Options: 'AMBER99SB_ILDN' | 'AMBER03' | 'CHARMM27' | 'OPLSAA' | 'GROMOS43a1'
                 | ...

            References
            ---------
            .. [5] Vitalini, F., Noe, F. and Keller, B. (2015):
                A basis set for peptides for the variational approach to conformational kinetics. (In review).

        """
        self.order = order
        self.ff = ff
        self.RBV = np.load('ResiduesEigenvectors/Ac_'+restype+'_NHMe.npz')
        self.RBV = self.RBV[ff][:,:order+1]
       
    def map(self, X):
        """
        Parameters
        ---------
        X : ndarray(T, 2)
            phi / psi angles for one residue in a  trajectory of length T
	    The array is angles in degrees between 0 and 360

        Returns
        -------
        Y : ndarray(T, order+1)
            The basis set values, evaluated for every trajectory time point and every function order.
            The first column is always the trivial 1-function.

        """
	binwidth=10 #binwidth
	i=np.floor(X[:,0]/binwidth) #binned phi angle
	j=np.floor(X[:,1]/binwidth) #binned psi angle
	n=i*(nbins)+j #n_microstates traj
        ev_traj=self.RBV[n,:] #project n_microstate trajectory of the monomer onto the eigenvector as ev_traj[t]=RBV[n[t]]
       	return ev_traj


class RamachandranProductBasis(BasisSet):

    def __init__(self, sequence, include_res=None, max_order=2, n_excite=2, ff='AMBER99SB_ILDN'):
        """Basis set as vector product of the residue centered basis functions evaluated on the phi/psi torsion space of the single amino acids.

        Parameters
        ---------
        sequence : iterable of str
            Residue types along the chain as single letter code, e.g. ['A', 'G', 'A', 'G', 'V']
	include_res : boolian array, defaule=None
	    Residues to include along the chain, eg ['True', 'True', 'False', 'False', 'True']. 
	    If not specified all residues are used.
        max_order : int, default=2
            maximum order of basis function family. Here equal to the number of nontrivial MSM eigenvectors
            on phi/psi space to be used.
        n_excite : int, default=2
            number of Ramachandran basis functions allowed to be excited. The default value (2) will allow only single and double excitations as suggested in [5].
            Careful, all possible combinations (n_excite=len(sequence)) can be unmanageable. The minimal excitation is 1, where
            the excited states up to max_order are used for each amino acid in the chain.
        ff : str, default='AMBER99SB_ILDN'
	    Force field to be used. Options: 'AMBER99SB_ILDN' | 'AMBER03' | 'CHARMM27' | 'OPLSAA' | 'GROMOS43a1'

        References
        ---------
        .. [5] Vitalini, F., Noe, F. and Keller, B. (2015):
            A basis set for peptides for the variational approach to conformational kinetics. (In review).

        """

        self.sequence = sequence
	self.include_res = include_res
        self.max_order = max_order
        self.n_excite = n_excite
        self.ff = ff

#	ind = np.arange(len(self.sequence), dtype=int) #array of indices
#	self.sequence = [self.sequence[i] for i in ind[[self.include_res]] # new sequence including only the residues labled as true
	compress_ind = it.compress(self.sequence, self.include_res) # itertool object including only the residues labled as true
	self.sequence = [item for item in compress_ind] # new sequence including only the residues labled as true
	self.nres = len(self.sequence) # number of residues to include
	combinations = it.combinations(np.arange(self.nres), self.n_excite) # select residues to excite
	basis_set_list = [] # initialize list of basis functions
	basis_set_list.append(np.zeros(self.nres)) # add the basis function where none of the residues is excited
	for c in combinations: # loop over combinations
		current_basis = np.zeros(self.nres, dtype=int) # initialize basis array
		excitation_res = it.combinations_with_replacement(np.arange(self.max_order+1), self.n_excite) # all possible combinations of eigenvectors in the excited states
		counter=0 # set to avoid adding the 00000 basis function all the times
		for p in excitation_res: # loop over the combinations of eigenvectors
			current_basis[np.array(c)]=np.array(p) # select spots to excite
			if counter > 0:
				basis_set_list.append(current_basis) # add excitations to list
			counter = counter + 1
	self.num_prod=len(basis_set_list) # number of products
	self.basis_set_list=basis_set_list # make list of basis functions available out of the function

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
	X=X[:,np.repeat(self.include_res, 2)] # select only residues labled as true
	Prod=np.ones((X.shape[0],self.num_prod)) # initialize products matrix
	q=0 # counter for the residue
	for r in range(self.nres): # loop over residues
		RBV=RamachandranBasis(self.sequence[r], order=self.max_order, ff=self.ff) # load eigenvectord for the residue
		RBV_traj=RBV.map(X[:,q:q+2]) # project trajectory over eigenvectors
		Prod=Prod*RBV_traj[:,self.basis_set_list[:,r]] # vectorizing over all possible states of excitations of residue r and multiplying the projected trajectory to the full product
		q =+ 2

	return Prod
        
