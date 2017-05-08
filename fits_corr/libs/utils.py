# Python 3 Imports
from __future__ import print_function
from __future__ import division

# Python 2 Imports
from os import path, walk
import numpy as np
from scipy.linalg import eigh, svd
from glob import glob
import warnings
import fitsio
from progressbar import ProgressBar, AdaptiveETA, Percentage, AnimatedMarker

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Internal Functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
def _get_int_dtype(max_val):
    """
    Internal function to return the least bits required to hold values
    up to <max_val> to preserve memory.
    
    """
    dtypes = map(lambda x: 'int' + x, ['8', '16', '32', '64'])
    values = np.asarray(map(lambda x: np.iinfo(x).max, dtypes))
    return dtypes[np.argmin((values - max_val) < 0)]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Public Functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
def progressbar(message='Progress: '):
    """
    Public utility to display a progress bar over iterables.
    
    """
    return ProgressBar(
        widgets=[
            message,
            ' ',
            Percentage(), 
            ' [', 
            AnimatedMarker(markers='-\\||//-'),
            '] ', 
            AdaptiveETA()
        ]
    )

def get_raw_dirs(root):
    """
    Scans recursively down the <root> path looking for
    parent directories containing .fits files which can
    be globbed on.
    
    Input
    -----
    root      : The source directory path.
    
    Output
    ------
    Generator object with directory paths.
    
    """
    return (root for root, dirname, files in walk(root) if len(files) > 1)

def get_on_off_pairs(directory, show_progress=True):
    """
    Collects all ON and OFF pairs of .fits files in the
    given <directory> as a list of tuples.
    
    Input
    -----
    directory : The parent directory containing .fits files.
    
    Output
    ------
    List of tuples (ON_FILE_PATH, OFF_FILE_PATH).
    
    """
    pattern = '_OFF.fits'
    files = glob(path.join(directory, '*.fits'))
    iterator = progressbar('Collecting Pairs: ')(files) if show_progress else files
    bag = {}
    for w in iterator:
        on = w[:-9] + '.fits' if w.endswith(pattern) else w
        if on not in bag:
            bag[on] = on[:-5] + pattern
    return bag.items()

def open_fits(fits_fp):
    """
    Open a .fits file and fixes any unconfigured headers.
    
    Input
    -----
    fits_fp   : The file path to the .fits file.
    
    Output
    ------
    The FITS data object of class <astropy.io.fits.hdu.image.PrimaryHDU>
    whose header and data objects can be extracted via the dot notation.
    
    """
    return fitsio.read(fits_fp)

def remap_data(data, min_val=0, max_val=255, dtype='i2'):
    """
    Remaps a given data matrix or image into the specified range given
    by <min_val> and <max_val> as optionally casts the data to the given
    <dtype>.
    
    Inputs
    ------
    data      : The input data matrix.
    min_val   : Lower bound of mapped range.
    max_val   : Upper bound of mapped range.
    dtype     : A valid NumPy datatype e.g. int, float, i4, f4, etc.
    
    Output
    ------
    Remapped data matrix.
    
    """
    lower = data.min()
    upper = data.max()
    val_range = max_val - min_val
    remapped_data = (((data - lower)/(upper - lower)) * val_range) + min_val
    
    if dtype is 'int' or dtype.startswith('i'):
        return np.round(remapped_data).astype(dtype)
    else:
        return remapped_data.astype(dtype)

def whiten_data(data):
    """
    Whitens the given data matrix by centering and sphering the data.
    
    Input
    -----
    data      : The input data matrix.
    
    Output
    ------
    Whitened data matrix.
    
    """
    return (data-data.mean())/data.std()
    
def compute_entropy(probs):
    """
    Computes the entropy given <probs> as the probability distribution using
    the function:
    
        H(A) = sum[ -probability * log2(probability) ]
    
    Input
    -----
    probs     : N-dimensional probability distribution matrix or array.
    
    Output
    ------
    The entropy of the distribution.
    
    """
    valid_probs = probs[probs.nonzero()]
    return np.sum(-valid_probs * np.log2(valid_probs))

def build_rolling_index(N, d, radius=3):
    """
    Generate index grid for accessing each patch of a rolling window
    as a vector.
    
    Inputs
    ------
    N         : Size of the 1st dimension.
    d         : Size of the 2nd dimension.
    radius    : Rolling window radius.
    
    Outputs
    -------
    (Y, X) index pair each with dimensions (N-2r)(M-2r) x (2r+1)^2.
    
    """
    span = (2 * radius)
    width = span + 1
    d_range = np.arange(width)
    y_range = N - span
    x_range = d - span
    dtype = _get_int_dtype(max(y_range, x_range))
    
    Y = np.zeros((y_range, width), dtype=dtype)
    Y += d_range + np.arange(y_range).reshape(-1, 1)
    Y = np.repeat(Y, width, axis=1)
    Y = np.repeat(Y, x_range, axis=0)
    
    X = np.zeros((x_range, width), dtype=dtype)
    X += d_range + np.arange(x_range).reshape(-1, 1)
    X = np.tile(X, (y_range, width))
    
    return (Y, X)

def pca_feature_vector(data):
    """
    Extract the top component explained by PCA given the <data> matrix
    with options for regularized PCA via <reg>.
    
    Inputs
    ------
    data      : Input data matrix to featurize.
    
    Output
    ------
    1-D feature vector corresponding to the top component explained by
    PCA or regularized PCA.
    
    """
    X = data - data.mean(0)
    _, _, V = eigh_extract_component(X)
    
    return np.dot(X, V)

def extract_top_component(data, reg=False, eps=1e-16):
    """
    Deprecated Function.
    """
    N, D = data.shape
    
    if N > D:
        C = np.dot(data.T, data)/N
        d, U = eigh(C, check_finite=False, overwrite_a=True)
        fv = U[:, -1]/np.sqrt(d[-1] + eps) if reg else U[:, -1]
    else:
        U, s, V = svd_extract_component(data)
        fv = V/(s + eps) * np.sqrt(N) if reg else V
    
    return fv

def svd_extract_component(data):
    """
    Deprecated Function.
    """
    U, s, V = svd(data, full_matrices=False, check_finite=False)
    
    return (U[:, 0], s[0], V[0])

def eigh_extract_component(data, compute_uv=False):
    """
    # TODO
    """
    N, D = data.shape
    U, s, V = None, None, None
    
    if N > D:
        C = np.dot(data.T, data)/N
        s, V = eigh(C, check_finite=False, overwrite_a=True)
        s, V = s[-1], V[:, -1]
        if compute_uv:
            U = np.dot(data, V)/np.sqrt(s * N)
    else:
        C = np.dot(data, data.T)/D
        s, U = eigh(C, check_finite=False, overwrite_a=True)
        s, U = s[-1], U[:, -1]
        V = np.dot(data.T, U)/np.sqrt(s * D)
    
    return U, s, V

def outer_product(data_a, data_b):
    """
    # TODO
    """
    vec_a = data_a.reshape(-1, 1)
    vec_b = data_b.reshape(1, -1)
    
    return np.dot(vec_a, vec_b)
    

def pcc_estimator(data_a, data_b):
    """
    # TODO
    """
    data_cov = np.sum(data_a * data_b)
    data_mags = np.sqrt(np.square(data_a).sum() * np.square(data_b).sum())
    
    return (data_cov/data_mags).__float__()

def prob_hist_2d(data_a, data_b, bins=256):
    """
    # TODO
    """
    N = data_a.size
    bin_r = np.arange(bins+1)
    n_w = np.ones((N,))
    
    low_x, high_x = data_a.min(), data_a.max()
    low_y, high_y = data_b.min(), data_b.max()
    
    h_x = ((high_x - low_x)/bins)
    h_y = ((high_y - low_y)/bins)
    
    edges_x = (bin_r * h_x) + low_x
    edges_y = (bin_r * h_y) + low_y
    
    edges_x[-1] = high_x + 1
    edges_y[-1] = high_y + 1
    
    i = np.searchsorted(edges_x, data_a, side='right') - 1
    j = np.searchsorted(edges_y, data_b, side='right') - 1
    
    hist = np.bincount(i*bins + j, weights=n_w, minlength=bins*bins)
    
    return hist.reshape((bins, bins))/N

def mutual_info(data_a, data_b, bins=256, norm=True):
    """
    Compute the mutual information via the joint entropy of <data_a>
    and <data_b>.
    
    Inputs
    ------
    data_a    : First data pair member.
    data_b    : Second data pair member.
    bins      : Type (int). Defaults to 256.
    norm      : Return normalized value of MI.
    
    Output
    ------
    (Normalized) Mutual Information score.
    
    """
    assert data_a.shape == data_b.shape, 'Error: Dimension mismatch!'
    
    joint_probs = prob_hist_2d(data_a, data_b, bins=bins)
    
    marginal_a = joint_probs.sum(1)
    marginal_b = joint_probs.sum(0)
    
    entropy_a = compute_entropy(marginal_a)
    entropy_b = compute_entropy(marginal_b)
    joint_entropy = compute_entropy(joint_probs)
    
    mutual_i = entropy_a + entropy_b - joint_entropy
    retval = (mutual_i/np.sqrt(entropy_a * entropy_b)) if norm else mutual_i
    
    return retval.__float__()