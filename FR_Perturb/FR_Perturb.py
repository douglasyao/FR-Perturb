'''
Functions for running FR-Perturb
'''
import numpy as np
import pandas as pd
import scipy 
import spams
from sklearn.linear_model import Lasso, LassoCV
from FR_Perturb.utils import *

def factorize_recover(exp_mat, p_mat_pd, rank, lambda1, log, alpha=None, iter=50, n_cv=10):
    ''' 
    Normalize expression matrix and run factorize-recover algorithm.

    Parameters:
    - exp_mat: numpy array, cell x gene raw expression matrix
    - p_mat_pd: pandas DataFrame, cell x perturbation indicator matrix (colnames must contain perturbation names)
    - rank: int, Rank parameter 
    - lambda1: float, Sparse PCA parameter (during factorize step)
    - lambda2: float, LASSO parameter (during recover step)
    - iter: Number of iterations to run sparse PCA

    Returns: 
    - B: numpy array, perturbation x gene unscaled LFC effect size matrix
    - U_tilde: numpy array, cell x module projected activity matrix 
    - U: numpy array, perturbation x module activity matrix
    - W: numpy array, module x gene dictionary matrix
    - p_mat: numpy array, cell x perturbation indicator matrix
    '''

    log.info('Running factorize step...  ')
    # exp_mat = exp_mat[:,np.squeeze(np.array(exp_mat.sum(axis = 0))) != 0]
    keep_cells = np.array(p_mat_pd.sum(axis = 1) > 0)[:,0]
    p_mat = p_mat_pd[keep_cells, :]
    # p_mat_pd = np.asfortranarray(p_mat_pd[keep_cells]).astype(np.float32)
    W = spams.trainDL(exp_mat, K=rank, lambda1=lambda1, iter=iter, verbose=False)
    U_tilde = spams.lasso(exp_mat, D=W, lambda1=lambda1, verbose=False)
    U_tilde = U_tilde[:, keep_cells]
    U_tilde = np.asfortranarray(U_tilde.T.todense()).astype(np.float32)
    W = W.T

    log.info('Running recover step... ')

    if alpha is None:
        U = LassoCV(n_alphas=10, fit_intercept=False)
        alphas = []
        for i in range(n_cv):
            U.fit(p_mat, U_tilde[:,i])
            alphas.append(U.alpha_)
        alpha = np.mean(alphas)
        
    U = Lasso(alpha = alpha, fit_intercept=False)
    U.fit(p_mat, U_tilde)
    U = U.coef_.T
    B = U.dot(W)
    
    return B, U, U_tilde, W, p_mat, alpha


def factorize_recover_large_dataset(mat, p_mat_pd, pnames, rank, n_batches, control_perturbation_name, log, alpha=None, cov_mat=None, n_cv=10):
    '''
    Normalize expression matrix and perform factorize-recover for a large dataset that cannot be fully stored in memory in dense format. 

    Parameters: 
        - mat: m x n sparse log(TP10K+1) expression matrix stored in CSR format. m = num cells, n = num genes.
        - p_mat_pd: m x p sparse perturbation indicator matrix. m = num cells, p = num perturbations. 
        - rank: number of principal components to compute
        - alpha: LASSO tuning parameter
        - n_batches: int, number of batches. At any given point, 1/n_batches of the # of entries in mat will be stored on disk in dense format.
        - control_perturbation_name: str, comma-delimited list of control perturbation names
        - cov_mat: (Optional) m x c pandas DataFrame of covariates. m = num cells, c = num covariates.

    Returns: 
        - B: numpy array, perturbation x gene unscaled LFC effect size matrix
        - U_tilde: numpy array, cell x module projected activity matrix 
        - U: numpy array, perturbation x module activity matrix
        - W: numpy array, module x gene dictionary matrix    
        - p_mat: numpy array, cell x perturbation indicator matrix
    '''
    
    if cov_mat is not None: 
        cov_mat = pd.get_dummies(cov_mat, drop_first=True) # Convert factors to dummy variables
        cov_means = cov_mat.values.mean(axis = 0) 
        cov_mat = cov_mat.values - cov_means[np.newaxis, :] # Center covariates
        cov_mat = np.c_[np.ones((cov_mat.shape[0], 1)), cov_mat] # Append intercept

    log.info('Getting mean control expression...   ')
    n_guides = p_mat_pd.sum(axis = 1)
    nt_names = control_perturbation_name.split(',')
    nt_idx = np.where(np.isin(pnames, nt_names))[0]
    
    if len(nt_idx) != len(nt_names):
        raise ValueError('Provided control perturbation names not found in data.')
    
    ctrl_idx = np.where(np.logical_and(n_guides == 1, p_mat_pd[:,nt_idx].sum(axis = 1) != 0))[0]

    if len(ctrl_idx) < 100:
        log.info('WARNING: Too few (<100) control cells, effect sizes will be computed relative to mean expression across all cells instead.')
        ctrl_idx = np.array(range(mat.shape[0]))
    if cov_mat is not None:
        ctrl_exp = regress_covariates_and_get_mean_expression_large_dataset(mat, cov_mat, n_batches, ctrl_idx)
    else:
        ctrl_exp = dat.X[ctrl_idx,:].mean(axis=0)

    keep_cells = np.array(p_mat_pd.sum(axis = 1) > 0)[:,0]
    p_mat = p_mat_pd[keep_cells, :]

    # factorize
    log.info('Running factorize step...  ')
    U_tilde, s, W = _normalize_and_svd_large_dataset(mat, rank, ctrl_exp, n_batches, cov_mat=cov_mat, num_iterations=3, log=log)
    U_tilde = U_tilde[keep_cells, :] * s
    U_tilde = np.asfortranarray(U_tilde).astype(np.float32)
    W = W.T

    # recover
    log.info('Running recover step...  ')
    if alpha is None:
        U = LassoCV(n_alphas=10, fit_intercept=False)
        alphas = []
        for i in range(n_cv):
            U.fit(p_mat, U_tilde[:,i])
            alphas.append(U.alpha_)
        alpha = np.mean(alphas)
    U = Lasso(alpha = alpha, fit_intercept=False)
    U.fit(p_mat, U_tilde)
    U = U.coef_.T
    B = U.dot(W)

    return B, U, U_tilde, W, p_mat, alpha


def _normalize_and_svd_large_dataset(mat, rank, ctrl_exp, n_batches, log, cov_mat=None, num_iterations=3):
    '''
    Normalize expression matrix and perform randomized SVD in batches as described in Halko et al. 
    "AN ALGORITHM FOR THE PRINCIPAL COMPONENT ANALYSIS OF LARGE DATA SETS".

    Parameters: 
        - mat: m x n sparse log(TP10K+1) expression matrix stored in CSR format. m = num cells, n = num genes.
        - rank: number of principal components to compute
        - ctrl_exp: numpy vector of length n. Mean normalized expression values in control cells
        - n_batches: int, number of batches. At any given point, 1/n_batches of the # of entries in mat will be stored on disk in dense format.
        - cov_mat: (Optional) m x c numpy matrix of centered covariates w/intercept. m = num cells, c = num covariates
        - num_iterations: int, number of iterations to run algorithm. 3 performs well in testing. 

    Returns: 
        - U, s, V: Eigenvectors and singular values
    '''

    # Normalize and factorize at the same time
    G = np.random.randn(mat.shape[1], rank + 2)
    log.info('Performing matrix multiplication 1 out of {}...   '.format(2*num_iterations + 2))
    H0 = _normalize_and_matmul_in_batches(mat=mat, mat2=G, n_batches=n_batches, cov_mat=cov_mat, ctrl_exp=ctrl_exp, log=log)
    prev_H = H0
    
    for i in range(num_iterations):
        log.info('Performing matrix multiplication {} out of {}...   '.format(2*i+2, 2*num_iterations + 2))
        AtHi_1 = _normalize_and_matmul_in_batches_transpose(mat=mat, mat2=prev_H, n_batches=n_batches, cov_mat=cov_mat, ctrl_exp=ctrl_exp, log=log)
        log.info('Performing matrix multiplication {} out of {}...   '.format(2*i+3, 2*num_iterations + 2))
        Hi = _normalize_and_matmul_in_batches(mat=mat, mat2=AtHi_1, n_batches=n_batches, cov_mat=cov_mat, ctrl_exp=ctrl_exp, log=log)
        prev_H = Hi
        H0 = np.c_[H0, Hi]

    Q, R, _ = scipy.linalg.qr(H0, mode='economic', pivoting=True)
    log.info('Performing matrix multiplication {} out of {}...   '.format(2*num_iterations + 2, 2*num_iterations + 2))
    T = _normalize_and_matmul_in_batches_transpose(mat=mat, mat2=Q, n_batches=n_batches, cov_mat=cov_mat, ctrl_exp=ctrl_exp, log=log)
    V_tilde, s_tilde, Wt = scipy.linalg.svd(T, full_matrices=False)
    U_tilde = Q.dot(Wt.T)
    U = U_tilde[:,:rank]
    V = V_tilde[:,:rank]
    s = s_tilde[:rank]

    # with open('temp.pickle', 'wb') as file:
    #     pickle.dump([U, s, V], file)

    return U, s, V

def _normalize_and_matmul_in_batches(mat, mat2, n_batches, cov_mat, ctrl_exp, log):
    '''
    Normalizes (i.e. regresses out covariates and centers on control expression) a sparse expression matrix on-the-fly, then multiplies the normalized matrix by 
    another matrix without ever storing the full dense matrix on disk. 
    
    Parameters:
        - mat: m x n sparse matrix expression stored in CSR format. m = num cells, n = num genes. 
        - mat2: n x k numpy matrix to multiply with mat. 
        - cov_mat: m x c numpy matrix of covariates to regress out of mat. m = num cells, c = num covariates. 
        - n_batches: int, number of batches. At any given point, 1/n_batches of the # of entries in mat will be stored on disk in dense format. 
        - ctrl_exp: numpy vector of length n. Mean normalized expression values in control cells
        - verbose: bool, whether to output progress

    Returns:
        - numpy array, m x k matrix 
    '''
    
    final_mat = np.zeros((mat.shape[0], mat2.shape[1]))
    batches = list(np.linspace(0, mat.shape[1], n_batches+1))
    batches = [int(x) for x in batches]
    
    for i in range(n_batches):
        start = batches[i]
        end = batches[i+1]
        log.info('Batch {} of {}...   '. format(i+1, n_batches))
        temp_mat = mat[:,start:end]

        if scipy.sparse.issparse(temp_mat):
            temp_mat = temp_mat.todense()
        
        if cov_mat is not None:
            temp_reg = scipy.linalg.lstsq(cov_mat, temp_mat, lapack_driver='gelsy')[0]
            temp_mat = temp_mat - cov_mat.dot(temp_reg)
            
        temp_mat = temp_mat - ctrl_exp[start:end]
        temp_mat = temp_mat.dot(mat2[start:end])
        final_mat += temp_mat
    return final_mat

def _normalize_and_matmul_in_batches_transpose(mat, mat2, n_batches, cov_mat, ctrl_exp, log):
    '''
    Same as normalize_and_matmul_in_batches(), but instead mat is transposed. 
    '''
    
    final_mat = np.empty((mat.shape[1], mat2.shape[1]))
    batches = list(np.linspace(0, mat.shape[1], n_batches+1))
    batches = [int(x) for x in batches]

    for i in range(n_batches):
        start = batches[i]
        end = batches[i+1]
        
        log.info('Batch {} of {}...   '. format(i+1, n_batches))        
        temp_mat = mat[:,start:end]

        if scipy.sparse.issparse(temp_mat):
            temp_mat = temp_mat.todense()
        
        if cov_mat is not None:
            temp_reg = scipy.linalg.lstsq(cov_mat, temp_mat, lapack_driver='gelsy')[0]
            temp_mat = temp_mat - cov_mat.dot(temp_reg)
        
        temp_mat = temp_mat - ctrl_exp[start:end]
        temp_mat = temp_mat.T.dot(mat2)
        final_mat[start:end] = temp_mat
    return final_mat