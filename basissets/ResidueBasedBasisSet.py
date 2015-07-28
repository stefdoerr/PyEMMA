__author__ = 'vitalini'


def ResidueBasedBasisSet(list_of_trajectories, prefix, ff='ff_AMBER99SB-ILDN', seq, included_res, order=2):
	""" Computes eigenvector basis functions as described in [5]_

	Parameters
	---------
	list_of_trajectories: 	list of npy-files.
				Each npy-file is a feature trajectory, as np-array, of shape (T,N), where T is the number of timesteps and N = 2 x sequence, i.e. phi / psi angles for each residue in the sequence

	prefix:	string.
		Commmon prefix to be used for all files to be produced.

	ff:	string, default='ff_AMBER99SB-ILDN'
	        Force field to be used. Options: 'ff_AMBER99SB-ILDN' | 'ff_AMBER03' | 'ff_CHARMM27' | 'ff_OPLS-AA' | 'ff_GROMOS43a1'

	seq:	string.
		Peptide sequence in one-letter code. Example: deca-alanine seq=AAAAAAAAAA

	included_res:	string.
			String of 0 and 1 used to select which residues of the sequence to include in the model. 0 = residue not included; 1=residue included. Example: deca-alanine, model exluding even-numbered residues 1010101010

	order:	int, default=2; minimum value=1  maximum value=length(seq)
		Number of excitations to be included.
		If 1 only products up to one residue not at equilibrium are included. Example generic trimer: list of products {000 001 002 010 020 100 200 }
		If 2, only products up to two residues not at equilibrium are included. Example generic trimer: list of products {000 001 002 010 020 100 200 011 012 021 022 110 120 210 220 101 102 201 202}
		if 3, only products up to three residues not at equilibrium are included. Example generic trimer: list of products {000 001 002 010 020 100 200 011 012 021 022 110 120 210 220 101 102 201 202 111 112 121 122 211 212 221 222} 
	

	Returns
	-------
	Y:	List of lists of filenames, where the evaluation of all basis functions can be found. The files will be called prefix_trajnum_funcnum.npy, where trajnum is the trajectory number and funcnum is the number associated to the basis function, here a string indicating the residues products.
		Each of these files contains the evaluation of the function funcnum over the trajectory trajnum. They are saved by the function as np-array of shape (T,1), with T number of timesteps.

	References
	---------
	.. [5] Vitalini, F., Noe, F. and Keller, B. (2015):
	A basis set for peptides for the variational approach to conformational kinetics. (In review).

	"""
	# TODO:
	# identify list of all possible products
	# access table data corresponding to the chosen ff
	# read function values from the discrete phi/psi boxes corresponding to X (careful: must be vectorized over time - otherwise slow)
	# project trajectory over residue-based functions.
	# evaluate products
	# Save projected trajectory over basis functions
	# return Y
    pass

	# VARIABLES
	ResidueBasisSetPath='/group/ag_cmb/scratch/vitalini/BasisSetLibrary/'
