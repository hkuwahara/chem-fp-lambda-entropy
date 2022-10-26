import numpy as np
from sklearn.base import TransformerMixin


class EigenEntropyThreshold(TransformerMixin):
    """Impementation of the "Eigenvalue based entropy" method for comparing fontributions of differentfingerprint bits.
    Inherets from `sklearn.base.TransformerMixin` to allow easy integration into `sklearn` pipelines.

    Attributes
    ----------
    threshold : float
        Reduced level theshold parameter for determining a bit's contribution.
        Specified at init but can be altered when calling `transform`.
        Must be valued [0, 1)

    jitter : float
        Noise term to prevent division by zero errors when taking calculating the logarithm.

    deltas : np.ndarray[float]
        Difference in eigenvalue entropy between reference value (all bits) and specified bit removed.

    mask : np.ndarray[bool]
        Boolean mask indicating if a particular fingerprint bit meets the Eigenvalue entropy threhold or not.

    Notes
    -----
    The algorithm is implemented as described in reference [1] with one notable exception as the original paper
    performed `n` dot products (i.e. 1 for each bit in the fingerprint being assessed).
    Here a single dot product is calculated for the entire feature matrix when `fit` is called, afterwhich this
    result is sliced to save having to perform multiple dot product calculations, hence allowing it to scale better.

    References
    ----------
    [1] J Cheminform 13, 27 (2021). https://doi.org/10.1186/s13321-021-00506-2
    [2] https://github.com/hkuwahara/chem-fp-lambda-entropy  (no license provided, original implementation in `R`)
    """

    def __init__(self, threshold=0.1, jitter=1e-10, mask=None):
        """
        Parameters
        ----------
        threshold : float
            default = 0.1
            Reduced level theshold parameter for determining a bit's contribution.
            Specified at init but can be altered when calling `transform`.
            Must be valued [0, 1)

        jitter : float
            default = 1e-10 (value used in original paper and hence highly reccomended(
            Noise term to prevent division by zero errors when taking calculating the logarithm.

        mask : np.ndarray[bool]
            default = None
            Pre-computed mask from prior transformation, allows for masks to be re-used and save long, repeated
            computations for large feature matrices.

        Notes
        -----
        If a pre-computed mask is passed, then passing a new threshold to `transform` will raise a `ValueError` as the
        relevant computations will not have been performed when using a pre-computed array.

        """
        self.threshold = self._ensure_valid_threshold(threshold)
        self.jitter = abs(float(jitter))
        self.deltas = np.array([])

        if mask is None:
            self.base_mask = np.array([])
            self._is_fitted = False
            self._is_precomputed = False
        else:
            self.base_mask = self._ensure_valid_mask(mask)
            self._is_fitted = True
            self._is_precomputed = True

        self.mask = np.array([])

    @staticmethod
    def _ensure_valid_threshold(threshold):
        """Ensures that threshold value is float valued valued [0, 1).
        """
        threshold = abs(float(threshold))
        assert 0 <= threshold < 1.0, F'Threshold must be 0 <= threshold < 1, not {threshold}'
        return threshold

    @staticmethod
    def _ensure_valid_mask(mask):
        """Ensures that mask is a boolean row vector.
        """
        try:
            mask = np.asarray(mask).astype(bool)
            assert mask.ndim == 1
        except (ValueError, AssertionError):
            raise ValueError('Mask must be boolean row vector.')
        return mask

    @staticmethod
    def _calc_eigen_entropy(X, jitter):
        """Calculates the Eigenvalue based entropy for a given matrix `X`.

        Parameters
        ----------
        X : np.ndarray
            shape(n_entries, n_features)
            Feature matrix.

        jitter : float
            Noise term to prevent division by zero errors when taking calculating the logarithm.

        Returns
        -------
        entropy : float
            Eigenvalue based entropy for passed feature matrix.
        """
        eig = np.linalg.svd(X, compute_uv=False) ** 2
        normed_eig = eig / eig.sum()
        entropy = -np.sum(normed_eig * np.log10(normed_eig + jitter)) / np.log10(X.shape[1])
        return entropy

    @staticmethod
    def _determine_mask(deltas, threshold):
        """Deteremines which bits meet the condition for inclusion based on entropy calculations.

        Parameters
        ----------
        deltas : np.ndarray[float]
            Difference between reference entropy and determined entropies in absence of each bit

        threshold : float
            Reduced level theshold parameter for determining a bit's contribution.
            Specified at init but can be altered when calling `transform`.
            Must be valued [0, 1)

        Returns
        -------
        mask : np.ndarray[bool]
            shape(len(deltas), )
            Boolean mask indicating if a particular fingerprint bit meets the Eigenvalue entropy threhold or not.
        """
        mask = np.abs(deltas) >= (threshold * np.sqrt(np.mean(deltas ** 2)))
        return mask

    def _determine_eigen_spectrum(self, A):
        """Determine the eigenvalue based entropy for each column in `A`.

        Parameters
        ----------
        A : np.ndarray
            shape(n_entries, n_features)
            Feature matrix.

        Returns
        -------
        hi : np.ndarray[float]
            Eigen based entropy values for each feature in `A`.
        """
        n = A.shape[1]
        n_bits = np.arange(n)
        hi = np.zeros(n)

        for bit in n_bits:
            bits_assessed = np.delete(n_bits, bit)
            ai = A[:, bits_assessed]
            hi[bit] = self._calc_eigen_entropy(ai, self.jitter)

        return hi

    def fit(self, X):
        """Determine the feature mask for the passed featuer matrix based on the threshold set at initialisation.

        Parameters
        ----------
        X : np.ndarray
            shape(n_entries, n_features)
            Feature matrix.

        Returns
        -------
        None : updates {self.deltas, self.mask}
        """
        A = np.asarray(X)
        h0 = self._calc_eigen_entropy(A, self.jitter)
        hi = self._determine_eigen_spectrum(A)

        self.deltas = hi - h0
        self.base_mask = self._determine_mask(self.deltas, self.threshold)
        self._is_fitted = True
        return self

    def transform(self, X, threshold=None):
        """Remove columns from passed feature matrix `X` based on determined `mask`.

        Notes
        -----
        A different threshold can be passed here which will determine a "new" mask to be used during the function
        call only, i.e. the `mask` attribute will not be updated when passing a threshold value.
        This is to allow for the evaluation of the mask at different threshold values without recalculating the eigen
        entropy.

        Notes
        -----
        Calling this method sets the `mask` attribute which allows for examination of which features are being
        masked or not.
        These are set as an attribute rather than returned when calling `transform` to maintain consistency with the
        `sklearn` transformer API.

        Raises
        ------
        ValueError : if using a pre-computed mask but attempt to specify a new threhold when calling `transform`.

        Parameters
        ----------
        X : np.ndarray
            shape(n_entries, n_features)
            Feature matrix.

        threshold : float
            Reduced level theshold parameter for determining a bit's contribution.
            Specified at init but can be altered when calling `transform`.
            Must be valued [0, 1)

        Returns
        -------
        X_masked : np.ndarray
            shape(n_entries, n_unmasked_features)
            Feature matrix.
        """
        assert self._is_fitted, 'Fit method must be called or precomputed mask specified before transforming `X`.'

        if threshold is not None:
            if self._is_precomputed:
                raise ValueError('New thresholds cannot be specified when using pre-computed feature masks.')
            else:
                threshold = self._ensure_valid_threshold(threshold)
                self.mask = self._determine_mask(self.deltas, threshold)
        else:
            self.mask = self.base_mask

        X = np.asarray(X)

        if len(self.mask) != X.shape[1]:
            raise ValueError(F'Mask of shape {self.mask.shape} has mismatched number of entries for `X` {X.shape[1]}.')

        return X[:, self.mask]
