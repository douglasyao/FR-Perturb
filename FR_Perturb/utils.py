'''
Other functions
'''
import pandas as pd
import numpy as np
import scipy
import statsmodels.api as sma
import statsmodels.stats as sms
import functools
    
def fit_skew_norm(t_star, t_nulls, side='both'):
    '''
    Compute p-values by fitting skew normal to null distribution. 
    
    Parameters
    ----------
    t_star: Test statistic
    t_null: Null statistics
    side: Which side to compute pvalues (left, right, or both)
    
    Returns
    ---------
    P-value
    '''
    
    if t_star == 0:
        p = 1
    else:
        fit = scipy.stats.skewnorm.fit(t_nulls)
        if side == 'left':
            p = scipy.stats.skewnorm.cdf(t_star, *fit)
        elif side == 'right':
            p = 1 - scipy.stats.skewnorm.cdf(t_star, *fit)
        elif side == 'both':
            p = scipy.stats.skewnorm.cdf(t_star, *fit)
            p = 2 * np.minimum(p, 1 - p)
        else:
            raise ValueError('Wrong side')
    return p

def scale_effs(B, logmeanexp, downsample_num = 25000, log_exp_baseline = 2):
    '''
    Scale effect sizes to mean expression using LOWESS. 
    
    Parameters
    ----------
    B: Perturbation x gene unscaled effect size matrix
    logmeanexp: Vector of log mean expression values to scale to
    downsample_num: Number of effects used to fit curve
    log_exp_baseline: Mean effect magnitude from this log expression is taken as the value to scale to
    
    Returns
    ---------
    B: Perturbation x gene scaled effect size matrix
    scale_factors: Per-gene scale factors 
    '''
    
    data_frac = min(1, downsample_num / np.prod(B.shape))
    
    if B.shape[1] != len(logmeanexp):
        raise ValueError('Number of genes differs')
    rand_idx = np.c_[np.random.randint(0, B.shape[0], downsample_num), 
                     np.random.randint(0, B.shape[1], downsample_num)]
    to_plot = np.c_[logmeanexp[rand_idx[:,1]], np.abs(B[rand_idx[:,0],rand_idx[:,1]])]
    to_plot = to_plot[np.where(to_plot[:,1] != 0)[0],:]
    to_plot[:,1] = np.log(to_plot[:,1])
    fit = sma.nonparametric.lowess(to_plot[:,1], to_plot[:,0], return_sorted=False, xvals = logmeanexp)
    baseline = fit[min(i for i,x in enumerate(logmeanexp) if x > log_exp_baseline)]
    scale_factors = np.exp(fit - baseline)
    B = B / scale_factors
    return B, scale_factors

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = functools.reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

def signif(X, n):
    '''Round elements of a pandas DF X to n significant figures'''
    def func(x):
        if x == 0:
            return 0
        else:
            return round(x, n - 1 - int(np.floor(np.log10(abs(x)))))
    return X.applymap(func) 

def convert_perturbation_vector_to_matrix(input, delimiter):
    '''
    Converts a pandas Series with multiple categories separated by a delimiter in each entry
    to a binary indicator matrix.

    Parameters:
    - categorical_data: pandas Series with string entries containing multiple categories separated by a delimiter.
    - delimiter: string, the delimiter used to separate multiple categories within each entry.

    Returns:
    - DataFrame: A binary indicator matrix where rows correspond to the index of the input series,
                 and columns correspond to the unique individual categories across all entries.
                 Each element in the matrix is a binary indicator (0 or 1), indicating the presence
                 of a category in each row.
    '''    
    split_categories = input.str.split(delimiter).apply(pd.Series).stack().reset_index(level=1, drop=True)
    indicator_matrix = pd.get_dummies(split_categories).groupby(level=0).max()
    return indicator_matrix

def regress_covariates(exp_mat, cov_mat):
    '''
    Regress out covariates from expression matrix. 
    Two orders of magnitude faster than scanpy.pp.regress_out function.
    
    Parameters:
        - exp_mat: numpy array, cell x gene raw expression matrix
        - cov_mat: pandas DataFrame, cell x covariate matrix

    Returns:
        - numpy array, cell x gene residual expression matrix
    '''
    
    cov_mat = pd.get_dummies(cov_mat, drop_first=True) # Convert factors to dummy variables
    cov_means = cov_mat.values.mean(axis = 0) 
    cov_mat = cov_mat.values - cov_means[np.newaxis, :] # Center covariates
    cov_mat = np.c_[np.ones((cov_mat.shape[0], 1)), cov_mat] # Append intercept
    
    if scipy.sparse.issparse(exp_mat):
        exp_mat = exp_mat.todense()

    lmfit = scipy.linalg.lstsq(cov_mat, exp_mat, lapack_driver='gelsy')[0]
    resids = exp_mat - cov_mat.dot(lmfit)
    return resids

def regress_covariates_and_get_mean_expression_large_dataset(mat, cov_mat, n_batches, idxs):
    '''
    Regresses out covariates from each column (gene) of an expression matrix, while only storing the entries for a subset of rows (cells). 
    This function is used when the original sparse expression matrix is extremely large and cannot be fully stored in memory in dense form.

    Covariates are regressed out one gene (column) at a time over all cells (rows) simultaneously. This function will do this for 
    batches of columns while only storing the values for a subset of rows.
    
    Parameters:
        - mat: m x n sparse matrix stored in CSR format
        - cov_mat: m x c numpy matrix of centered covariates w/intercept. m = num cells, c = num covariates
        - n_batches: int, number of column batches to perform regression over
        - idxs: 1-D numpy array or list, row indices of mat to retain

    Returns:
        - numpy array, len(idx)
    '''
    
    col_batch_size = round(mat.shape[1] / n_batches)
    final_mat = np.empty(mat.shape[1])
    for i in range(0, mat.shape[1], col_batch_size):
        start = i
        end = min(mat.shape[1], i+col_batch_size)
        temp_mat = mat[:,start:end]

        if scipy.sparse.issparse(temp_mat):
            temp_mat = np.array(temp_mat.todense())
            
        temp_reg = scipy.linalg.lstsq(cov_mat, temp_mat, lapack_driver='gelsy')[0]
        temp_mat = temp_mat - cov_mat.dot(temp_reg)
        temp_mat = temp_mat[idxs]
        final_mat[start:end] = temp_mat.mean(axis=0)
    return final_mat

def compute_correlation(mat1, mat2, idxs, sort, cutoffs):
    '''
    Compute correlation between two matrices
    '''
    mat1 = mat1.flatten()
    mat2 = mat2.flatten()
    sort = sort.flatten()
    sort_idxs = np.argsort(-np.abs(sort))

    cutoffs = np.round(np.array(cutoffs) * len(sort)).astype(int)
    corrs = []
    sign_concords = []
    for cutoff in cutoffs:
        corrs.append(np.corrcoef(mat1[sort_idxs[:cutoff]], mat2[sort_idxs[:cutoff]])[0,1])
        sign_concords.append(sign_concord(mat1[sort_idxs[:cutoff]], mat2[sort_idxs[:cutoff]]))
    return corrs, sign_concords

def sign_concord(x, y):
    '''
    Compute sign concordance of two vectors x and y. 0s in either x or y will have 0.5 sign concordance by default. 
    '''
    all_dat = pd.DataFrame({'x': x, 'y': y})
    all_dat.dropna(inplace=True)
    zeros = all_dat[(all_dat['x'] == 0) | (all_dat['y'] == 0)].index    
    nonzero_dat = all_dat.drop(zeros)
    nonzero_sc = np.mean(np.sign(nonzero_dat['x']) == np.sign(nonzero_dat['y']))
    if len(all_dat) > 0:
        sign_concordance = (nonzero_sc * len(nonzero_dat) + 0.5 * len(zeros)) / len(all_dat)
    else:
        sign_concordance = 0  # handle case with no valid data
    
    return sign_concordance

