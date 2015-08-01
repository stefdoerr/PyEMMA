from base import BasisSet
import numpy as np
import itertools as it
import os


class RamachandranBasis(BasisSet):

    def __init__(self, restype, order=2, ff='ff_AMBER99SB_ILDN', radiants=True):
        """
        Basis set on the phi/psi torsion space of a single amino acid.

        Parameters
        ---------
        restype :   str
                    Residue type, single letter code, e.g. alanine = 'A'.
        order : int, default=2
                order of basis function family. 
                Here equal to the number of nontrivial MSM eigenvectors on phi/psi space
                to be used. The output will have order+1 columns
        ff : str, default='ff_AMBER99SB_ILDN'
             Force field to be used. Options: 'ff_AMBER99SB_ILDN' | 'ff_AMBER03' | 'ff_CHARMM27' | 'ff_OPLSAA' | 'ff_GROMOS43a1'
             | ...
        radiants :  True
                    expected inputs as radiants. If false, treats input as degrees

        References
        ---------
        .. [5] Vitalini, F., Noe, F. and Keller, B. (2015):
        A basis set for peptides for the variational approach to conformational kinetics. (In review).
        """

        path_local = os.path.dirname(os.path.abspath(__file__))
        self.order = order
        self.ff = ff
        self.radiants = radiants
        self.RBV = np.load(path_local+'/ResiduesEigenvectors/Ac_'+restype+'_NHMe.npz')
        self.RBV = self.RBV[ff][:,:order+1]
       
    def map(self, X):
        """
        Parameters
        ---------
        X : ndarray(T, 2)
            phi / psi angles for one residue in a  trajectory of length T
            The array is angles in radiants between 0 and 360

        Returns
        -------
        Y : ndarray(T, order+1)
            The basis set values, evaluated for every trajectory time point and every function order.
            The first column is always the trivial 1-function.
        """
        binwidth=10 #binwidth
        nbins=360/binwidth # number of bins per coordinate
        if self.radiants==True:
            X=X*(180/np.pi) # convert in radiants
        if np.max(X)<360:
            X=X+180
        X[np.nonzero(X[:,0]>=360),0]=0
        X[np.nonzero(X[:,1]>=360),1]=0
        X[np.nonzero(X[:,0]<0),0]=0
        X[np.nonzero(X[:,1]<0),1]=0
        i=np.floor(X[:,0]/binwidth) #binned phi angle
        j=np.floor(X[:,1]/binwidth) #binned psi angle
        n=i*(nbins)+j #n_microstates traj
        n=np.array(n, dtype=int)
        ev_traj=self.RBV[n,:] #project n_microstate trajectory of the monomer onto the eigenvector as ev_traj[t]=RBV[n[t]]
        
        return ev_traj


class RamachandranProductBasis(BasisSet):

    def __init__(self, sequence, include_res=None, max_order=2, n_excite=2, ff='ff_AMBER99SB_ILDN', radiants=True):
        """
        Basis set as vector product of the residue centered basis functions evaluated 
        on the phi/psi torsion space of the single amino acids.

        Parameters
        ---------
        sequence :  iterable of str
                    Residue types along the chain as single letter code, 
                    e.g. ['A', 'G', 'A', 'G', 'V']
        include_res :   boolian array, defaule=None
                        Residues to include along the chain, 
                        eg ['True', 'True', 'False', 'False', 'True']. 
                        If not specified all residues are used.
        max_order :     int, default=2
                        maximum order of basis function family. 
                        Here equal to the number of nontrivial MSM eigenvectors
                        on phi/psi space to be used.
        n_excite :  int, default=2
                    number of Ramachandran basis functions allowed to be excited. 
                    The default value (2) will allow only single and double excitations 
                    as suggested in [5].
                    Careful, all possible combinations (n_excite=len(sequence)) can be unmanageable. 
                    The minimal excitation is 1, where the excited states up to max_order 
                    are used for each amino acid in the chain.
        ff :    str, default='ff_AMBER99SB_ILDN'
                Force field to be used. Options: 'ff_AMBER99SB_ILDN' | 'ff_AMBER03' | 'ff_CHARMM27' | 'ff_OPLSAA' | 'ff_GROMOS43a1'

        radiants :  True
                    expected inputs as radiants. If false, treats input as degrees


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
        self.radiants = radiants
	
        if self.include_res == None:
            self.include_res=np.ones(len(self.sequence), dtype=bool)

        compress_ind = it.compress(self.sequence, self.include_res) # itertool object including only the residues labled as true
        self.sequence = [item for item in compress_ind] # new sequence including only the residues labled as true
        self.nres = len(self.sequence) # number of residues to include
        combinations = it.combinations(np.arange(self.nres), self.n_excite) # select residues to excite
        basis_set_list = [] # initialize list of basis functions
        for c in combinations: # loop over combinations
            current_basis = np.zeros(self.nres, dtype=int) # initialize basis array
            excitation_res = it.product(np.arange(self.max_order+1), repeat=self.n_excite)# all possible combinations of eigenvectors in the excited states
            for p in excitation_res: # loop over the combinations of eigenvectors
                current_basis[np.array(c)]=np.array(p) # select spots to excite
                basis_set_list.append(np.copy(current_basis)) # add excitations to list
        basis_set_list=np.array(basis_set_list, dtype=int)
        #find unique elements in basis set list
        b = np.ascontiguousarray(basis_set_list).view(np.dtype((np.void, basis_set_list.dtype.itemsize * basis_set_list.shape[1])))
        _, idx =np.unique(b, return_index=True)
        unique_basis=basis_set_list[idx]
        self.num_prod=len(unique_basis) # number of products
        self.basis_set_list=np.array(unique_basis, dtype=int) # make list of basis functions available out of the function

    def map(self, X):
        """
        Parameters
        ---------
        X : ndarray(T, 2R)
            phi / psi angles for each of R residues for a trajectory of length T

        Returns
        -------
        Prod :  ndarray(T, n)
                The basis set values, evaluated for every trajectory time point.
                Each of the n elements represents a product of the residue basis
                sets.

        babasis_set_list :  nparray(nprod, nres)
                            list of basis functions as numpy array with 
                            nprod = nuber of basis functions and 
                            nres = number of residues
        """
        
        X=X[:,np.repeat(self.include_res, 2)] # select only residues labled as true
        Prod=np.ones((X.shape[0],self.num_prod)) # initialize products matrix
        q=0 # counter for the residue
        for r in range(self.nres): # loop over residues
            RBV=RamachandranBasis(self.sequence[r], order=self.max_order, ff=self.ff, radiants=self.radiants) # load eigenvectord for the residue
            RBV_traj=RBV.map(X[:,q:q+2]) # project trajectory over eigenvectors
            Prod=Prod*RBV_traj[:,self.basis_set_list[:,r]] # vectorizing over all possible states of excitations of residue r and multiplying the projected trajectory to the full product
            q = q + 2

        return Prod, self.basis_set_list
        
