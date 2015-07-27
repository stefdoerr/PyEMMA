__author__ = 'noe'

class BasisSet(object):

    def map(self, X):
        """
        Parameters
        ----------
        X : feature (d) or feature-trajectory (T x d)

        Result
        ------
        Y : basis (n) or (T x n)

        """
        raise NotImplementedError('must be implemented!')
