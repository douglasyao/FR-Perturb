#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy 
import spams
import scanpy
import time, sys, traceback, argparse
import os
import tqdm
import statsmodels as sm
import functools

class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module
    '''

    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        print(msg, file=self.log_fh)
        print(msg)
        
def regress_covariates(dat, cov_mat):
    '''
    Regress out covariates from expression matrix in place. 
    Two orders of magnitude faster than scanpy.pp.regress_out function.
    
    Parameters
    ----------
    dat: AnnData object
    cov_mat: Cell x covariate Dataframe
    '''
    
    cov_mat = pd.get_dummies(cov_mat, drop_first=True) # Convert factors to dummy variables
    cov_means = cov_mat.values.mean(axis = 0) 
    cov_mat = cov_mat.values - cov_means[np.newaxis, :] # Center covariates
    cov_mat = np.c_[np.ones((cov_mat.shape[0], 1)), cov_mat] # Append intercept
    
    if scipy.sparse.issparse(dat.X):
        dat.X = dat.X.todense()

    lmfit = scipy.linalg.lstsq(cov_mat, dat.X, lapack_driver='gelsy')[0]
    resids = dat.X - cov_mat.dot(lmfit)
    dat.X = resids
    
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

# from matplotlib import pyplot as plot
# downsample_num = 25000
# log_exp_baseline = 2
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
    to_plot = np.c_[logmeanexp[rand_idx[:,1]], np.log(np.abs(B[rand_idx[:,0],rand_idx[:,1]]))]
    to_plot = to_plot[~np.isinf(to_plot[:,1]),:]
    fit = sm.api.nonparametric.lowess(to_plot[:,1], to_plot[:,0], return_sorted=False, xvals = logmeanexp)
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

os.environ['KMP_WARNINGS'] = 'off'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* Factorize-Recover for Perturb-seq analysis (FR-Perturb)\n"
MASTHEAD += "* by Douglas Yao 2022 \n"
MASTHEAD += "*********************************************************************\n"

parser = argparse.ArgumentParser()

### Required flags
parser.add_argument('--input-h5ad', default=None, type=str,
                    help='TODO')
parser.add_argument('--input-perturbation-matrix', default=None, type=str,
                    help='TODO')
parser.add_argument('--control-perturbation-name', default=None, type=str,
                    help='Comma-separated list of perturbation names that represent control perturbations')
parser.add_argument('--out', default=None, type=str,
                    help="Output prefix (including directory) for effect sizes (and optionally p-values)")

### Optional
parser.add_argument('--compute-pval', default=False, action='store_true',
                    help='Whether or not to compute p-values for all effect size estimates by permutation testing')           
parser.add_argument('--rank', default=20, type=int,
                    help='Hyperparameter determining the rank of the matrix during the factorize step')
parser.add_argument('--lambda1', default=0.1, type=float,
                    help='Hyperparameter determining the sparsity of the factor matrix during the factorize step of the method. Higher value = more sparse.')
parser.add_argument('--lambda2', default=10, type=float,
                    help='Hyperparameter determining the sparsity of learned effects during the recover step of the method. Higher value = more sparse.')
parser.add_argument('--covariates', default=None, type=str,
                    help='Comma-separated list of covariate names to regress out of the expression matrix (names must match the column names in the meta-data of the Seurat object)')
parser.add_argument('--guide-overloaded', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from guide-overloading')
parser.add_argument('--droplet-overloaded', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from droplet-overloading')
parser.add_argument('--num-perms', default=10000, type=int,
                    help='Number of permutations when doing permutation testing')
parser.add_argument('--fit-zero-pval', default=False, action='store_true',
                    help='Compute p-values by fitting skew-t distribution to null distribution (allows for p-values below 1/num_perms, but significantly increases compute time)')
parser.add_argument('--output-factor-matrices', default=False, action='store_true',
                    help='Whether or not to output the latent gene expression factor matrices in addition to the full effect sizes matrix')

if __name__ == '__main__':

    args = parser.parse_args()
    log = Logger(args.out + '.log')
    
    defaults = vars(parser.parse_args(''))
    opts = vars(args)
    non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
    header = MASTHEAD
    header += "Call: \n"
    header += './run_FR_perturb.py \\\n'
    options = ['--' + x.replace('_', '-') + ' ' + str(opts[x]) + ' \\' for x in non_defaults]
    header += '\n'.join(options).replace('True', '').replace('False', '')
    header = header[0:-1] + '\n'
    log.log(header)
    log.log('Beginning analysis at {T}'.format(T=time.ctime()))
    start_time = time.time()
    
    if not args.input_h5ad:
        raise ValueError('Must specify --input-h5ad')
    if not args.input_perturbation_matrix:
        raise ValueError('Must specify --input-perturbation-matrix')
    if not args.control_perturbation_name:
        raise ValueError('Must specify --control-perturbation-name')
    if not args.out:
        raise ValueError('Must specify --out')
        
    if args.guide_overloaded and args.droplet_overloaded:
        raise ValueError('Only one of --guide-overloaded and --droplet-overloaded should be set')
    elif args.droplet_overloaded:
        overload_type = 'droplet'
    else:
        overload_type = 'guide'
        
    log.log('Loading input data...  ')
    dat = scanpy.read_h5ad(args.input_h5ad)
    p_mat_pd = pd.read_csv(args.input_perturbation_matrix, index_col = 0)
    if not dat.obs.index.equals(p_mat_pd.columns):
        raise ValueError('Cell names in perturbation matrix do not match cell names in expression matrix')
    log.log('Done')
    
    # center rows of p_mat
    if overload_type == 'droplet':
        guides_per_cell = np.array(p_mat_pd.sum(axis = 0))
        p_mat_pd = p_mat_pd.divide(guides_per_cell[np.newaxis, :])
        
    # get covariates
    if args.covariates:
        cov_names = args.covariates.split(',')
        cov_mat = dat.obs[cov_names]

    # regress out covariates 
    log.log('Regressing out covariates and centering expression matrix...  ')
    scanpy.pp.normalize_total(dat, target_sum = 10000)
    logmeanexp = np.log(np.array(dat.X.mean(axis = 0))[0,:])
    scanpy.pp.log1p(dat)
    if args.covariates:
        regress_covariates(dat, cov_mat)
    else:
        dat.X = dat.X - dat.X.mean(axis = 0)
    
    # center expression matrix based on control expression
    n_guides = p_mat_pd.values.sum(axis = 0)
    nt_names = args.control_perturbation_name.split(',')
    ctrl_idx = np.logical_and(n_guides == 1, p_mat_pd.loc[nt_names].sum(axis = 0).values != 0)
    ctrl_exp = dat.X[ctrl_idx,:].mean(axis = 0)
    dat.X = dat.X - ctrl_exp
    log.log('Done')
    
    # Factorize
    log.log('Factorizing expression matrix... ')
    dat.X = dat.X[:,np.array(dat.X.sum(axis = 0))[0,:] != 0]
    keep_cells = p_mat_pd.sum(axis = 0) > 0
    p_mat = np.asfortranarray(p_mat_pd.loc[:, keep_cells].T).astype(float)
    W = spams.trainDL(np.asfortranarray(dat.X.T), K=args.rank, lambda1=args.lambda1, iter=50, verbose=False)
    U_tilde = spams.lasso(np.asfortranarray(dat.X.T), D=W, lambda1=args.lambda1, verbose=False)
    U_tilde = U_tilde[:, keep_cells]
    U_tilde = np.asfortranarray(U_tilde.T.todense())
    W = W.T
    log.log('Done')
    
    # Recover
    log.log('Regressing left factor matrix on perturbation design matrix...  ')
    U = spams.lasso(U_tilde, D=p_mat, lambda1=args.lambda2, verbose=False)
    B = U.dot(W)
    log.log('Done')
    
    # Compute pvalues by permutation testing
    if args.compute_pval:
        if not args.fit_zero_pval:
            log.log('Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
            pvals = np.zeros((B.shape))
            for i in tqdm.tqdm(range(args.num_perms)):
                p_mat_perm = np.asfortranarray(p_mat[np.random.permutation(p_mat.shape[0]),:])
                U_perm = spams.lasso(U_tilde, D=p_mat_perm, lambda1=args.lambda2, verbose=False)
                B_perm = U_perm.dot(W)
                temp_indices = B < B_perm
                pvals[temp_indices] = pvals[temp_indices] + 1
            pvals /= args.num_perms
            pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
            pvals *= 2 
            pvals = (pvals * args.num_perms + 1) / (args.num_perms + 1)
            log.log('Done')
        else:
            args.num_perms = 500
            log.log('Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
            B_perms = np.empty((np.product(B.shape), args.num_perms))
            for i in tqdm.tqdm(range(args.num_perms)):
                p_mat_perm = np.asfortranarray(p_mat[np.random.permutation(p_mat.shape[0]),:])
                U_perm = spams.lasso(U_tilde, D=p_mat_perm, lambda1=args.lambda2, verbose=False)
                B_perm = U_perm.dot(W)
                B_perms[:,i] = np.ravel(B_perm)
            pvals = (B_perms < np.ravel(B)[:,np.newaxis]).sum(axis=1) / B_perms.shape[1]
            pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
            pvals *= 2 
            
            log.log('Fitting skew-normal distribution to effects with p=0 ({} total effects)...  '.format(np.sum(pvals == 0)))
            zero_indices = np.where(pvals == 0)[0]
            B_flattened = np.ravel(B)
            fitted_pvals = []
            for i in tqdm.tqdm(zero_indices):
                t_star = B_flattened[i]
                t_nulls = B_perms[i,:]
                if t_star == 0:
                    p = 1
                else:
                    p = fit_skew_norm(t_star, t_nulls)
                fitted_pvals.append(p)
            pvals[zero_indices] = fitted_pvals
            pvals = np.reshape(pvals, B.shape)
            log.log('Done')

        qvals = sm.stats.multitest.multipletests(pvals.flatten(), method = 'fdr_bh')[1]
        qvals = np.reshape(qvals, pvals.shape)
        pvals = pd.DataFrame(data = np.transpose(pvals), index = dat.var.features, columns = p_mat_pd.index)
        qvals = pd.DataFrame(data = np.transpose(qvals), index = dat.var.features, columns = p_mat_pd.index)
    
    log.log('Scaling effects...  ')
    B,_ = scale_effs(B, logmeanexp)
    log.log('Done')
    log.log('Outputting results...  ')
    B = pd.DataFrame(data = np.transpose(B), index = dat.var.features, columns = p_mat_pd.index)
    B.to_csv(args.out + '_LFCs.txt', sep = '\t')
    if args.compute_pval:
        pvals.to_csv(args.out + '_pvals.txt', sep = '\t')
        qvals.to_csv(args.out + '_qvals.txt', sep = '\t')
    log.log('All done!')       
    log.log('Analysis finished at {T}'.format(T=time.ctime()))
    time_elapsed = round(time.time() - start_time, 2)
    log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
        

