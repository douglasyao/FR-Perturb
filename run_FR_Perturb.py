#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy
import time, sys, traceback, argparse
import os
import tqdm
import random
import logging
import warnings
import pickle
from tqdm.contrib.concurrent import thread_map
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.linear_model import Lasso, LassoCV, ElasticNet
from FR_Perturb.FR_Perturb import *
from FR_Perturb.utils import *


# In[2]:


os.environ['KMP_WARNINGS'] = 'off'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* Factorize-Recover for Perturb-seq analysis (FR-Perturb)\n"
MASTHEAD += "* by Douglas Yao 2022 \n"
MASTHEAD += "*********************************************************************\n"

parser = argparse.ArgumentParser()

################ Required flags ##################
parser.add_argument('--input-h5ad', default=None, type=str,
                    help='h5ad file (from the AnnData package) containing raw gene expression counts for all cells')
parser.add_argument('--out', default=None, type=str,
                    help="Output prefix (including directory) for effect sizes")

# Either need to specify 
parser.add_argument('--input-perturbation-matrix', default=None, type=str,
                    help='Whitespace-delimited file containing a table with rows corresponding to cells and columns corresponding to perturbations. Cells containing a given perturbation should be indicated with a "1", otherwise "0".')
# Or
parser.add_argument('--perturbation-column-name', default=None, type=str,
                    help='Column name in the cell annotation matrix (stored in `.obs` in the h5ad object) that corresponds to which perturbations are in the cell.')
parser.add_argument('--perturbation-delimiter', default=None, type=str,
                    help='Delimiter used to separate perturbations in each cell (if the cell contains multiple perturbations). If the delimiter includes a double-dash --, then please enclose in both double and single quotes (e.g. '"--"'). This is a hack, sorry.')
###################################################


### Optional I/O
parser.add_argument('--covariates', default=None, type=str,
                    help='Comma-separated list of covariate names to regress out of the expression matrix (names must be found in the column names of `.obs` of the h5ad object)')
parser.add_argument('--control-perturbation-name', default=None, type=str,
                    help='If you want to compute effect sizes relative to mean expression in control cells rather than mean expression across all cells, specify a comma-separated list of control perturbation names here.')
parser.add_argument('--combine-control-perturbations', default=False, action='store_true',
                    help='Combine all separate control perturbations into a single "control" category when inferring effects.')
parser.add_argument('--output-factor-matrices', default=False, action='store_true',
                    help='Whether or not to output the latent gene expression factor matrices in addition to the full effect sizes matrix')
parser.add_argument('--gene-column-name', default=None, type=str,
                    help='Which column (by name) of `.var` in the h5ad object represents gene names. By default, the index is assumed to be gene names.')

# Interactions
parser.add_argument('--interaction-column-name', default=None, type=str,
                    help='Comma-separated list of variable names to compute interaction terms for (names must be found in the column names of `.obs` of the h5ad object)')

# Significance testing
parser.add_argument('--compute-pval', default=False, action='store_true',
                    help='Whether or not to compute p-values for all effect size estimates by permutation testing.')           
parser.add_argument('--num-perms', default=1000, type=int,
                    help='Number of permutations when doing permutation testing')
parser.add_argument('--fit-zero-pval', default=False, action='store_true',
                    help='Compute p-values by fitting skew-normal distribution to null distribution (allows for p-values below 1/num_perms, but significantly increases compute time)')
parser.add_argument('--multithreaded', default=False, action='store_true',
                    help='Use multithreading to fit skew-normal distributions.')
parser.add_argument('--temp-out', default=None, type=str,
                    help='Temporary file name, including directory (for batched permutation testing).')

# For running large datasets
parser.add_argument('--large-dataset', default=False, action='store_true',
                    help='Set this flag if dataset is large (>1M cells). Uses memory-efficient randomized partial SVD during factorize step of factorize-recover. --lambda1 parameter will be ignored.')
parser.add_argument('--batches', default=5, type=int,
                    help='Batch size when --large-dataset is set. Increasing this number reduces the amount of memory but increases running time.')

# Cross-validation                    
parser.add_argument('--cross-validate', default=False, action='store_true',
                    help='Perform 2-fold cross-validation on input data')
parser.add_argument('--cross-validate-runs', default=2, type=int,
                    help='Number of times to repeat cross-validation. Cross-validation accuracy is averaged across runs.')

# Hyperparameters
parser.add_argument('--rank', default=20, type=int,
                    help='Hyperparameter determining the rank of the matrix during the factorize step')
parser.add_argument('--spca-alpha', default=0.1, type=float,
                    help='Hyperparameter determining the sparsity of the factor matrix during the factorize step of the method. Higher value = more sparse.')
parser.add_argument('--lasso-alpha', default=0.0001, type=float,
                    help='Hyperparameter determining the sparsity of learned effects during the recover step of the method. Higher value = more sparse. If not specified, will determine automatically through cross-validation.')

# For guide/cell pooled experiments
parser.add_argument('--guide-pooled', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from guide pooling')
parser.add_argument('--cell-pooled', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from cell pooling')

# Inference type
parser.add_argument('--elastic-net', default=False, action='store_true',
                    help='Estimate effect sizes using elastic net (as done in Dixit et al. 2016 Cell) rather than FR-Perturb')
parser.add_argument('--elastic-net-genes-per-batch', default=None, type=int,
                    help='Split up the inference over batches of downstream genes to reduce memory usage. Input is the number of genes in each batch (larger = fewer batches but more memory)')


# In[ ]:


if __name__ == '__main__':

    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.FileHandler(args.out + ".log"),logging.StreamHandler()],
        level=logging.INFO,
    )
    log = logging.getLogger(__name__)

    try:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = MASTHEAD
        header += "Call: \n"
        header += './run_FR_perturb.py \\\n'
        options = ['--' + x.replace('_', '-') + ' ' + str(opts[x]) + ' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True', '').replace('False', '')
        header = header[0:-1] + '\n'
        log.info(header)
        log.info('Beginning analysis at {T}'.format(T=time.ctime()))
        start_time = time.time()

        if not args.input_h5ad:
            raise ValueError('Must specify --input-h5ad')
        if (args.input_perturbation_matrix is None) == (args.perturbation_column_name is None):
            raise ValueError('Must specify exactly one of --input-perturbation-matrix or --perturbation-column-name')
        if not args.out:
            raise ValueError('Must specify --out')
        if args.guide_pooled and args.cell_pooled:
            raise ValueError('Only one of --guide-pooled and --cell-pooled should be set')
        elif args.cell_pooled:
            overload_type = 'droplet'
        else:
            overload_type = 'guide' # use this version by default

        log.info('Loading input h5ad from {}...  '.format(args.input_h5ad))
        dat = scanpy.read_h5ad(args.input_h5ad)
        
        # get perturbation matrix
        if args.input_perturbation_matrix:
            log.info('Loading input perturbation matrix from {}...  '.format(args.input_perturbation_matrix))
            p_mat_pd = pd.read_csv(args.input_perturbation_matrix, index_col = 0, delim_whitespace=True)
            if not dat.obs.index.equals(p_mat_pd.index):
                raise ValueError('Cell names in perturbation matrix do not match cell names in expression matrix')
            pnames = p_mat_pd.columns.values
            p_mat_pd = scipy.sparse.csr_matrix(p_mat_pd.values)
        else:
            perts = dat.obs[args.perturbation_column_name]
            if isinstance(perts.dtype, pd.api.types.CategoricalDtype):
                if 'nan' not in perts.cat.categories: # replace empty perturbations with 'nan'
                    perts = perts.cat.add_categories(['nan'])
            perts[pd.isnull(perts)] = 'nan'
            mlb = MultiLabelBinarizer(sparse_output=True)

            if args.perturbation_delimiter is not None:
                if args.perturbation_delimiter[0] == '\"' and args.perturbation_delimiter[-1] == '\"': # if delimiter includes double dash
                    args.perturbation_delimiter = args.perturbation_delimiter[1:-1]
                p_mat_pd = mlb.fit_transform(perts.str.split(args.perturbation_delimiter, regex = False))
            else:
                p_mat_pd = mlb.fit_transform([[x] for x in perts.values])
            pnames=mlb.classes_
        
        # getting control perturbation info
        if args.control_perturbation_name is not None:
            pert_idx = np.where(np.isin(pnames, args.control_perturbation_name.split(',')))[0]
            if len(pert_idx) != len(args.control_perturbation_name.split(',')):
                raise ValueError('Provided control perturbation names not found in data. Please also double check that delimiter is specified correctly.')

            if args.combine_control_perturbations and len(args.control_perturbation_name.split(',')) > 1: # combine control perturbations into one category
                ctrl_perts = p_mat_pd[:,pert_idx]
                ctrl_perts = ctrl_perts.max(axis = 1)
                mask = np.ones(p_mat_pd.shape[1], dtype=bool)
                mask[pert_idx] = False
                p_mat_pd = scipy.sparse.hstack([p_mat_pd[:,mask], ctrl_perts], format='csr')
                pnames = np.append(pnames[mask], 'control')
                args.control_perturbation_name = 'control'

        # get covariates
        if args.covariates is not None:
            cov_names = args.covariates.split(',')
            cov_mat = dat.obs[cov_names]
            if np.sum(pd.isnull(cov_mat.values)) > 0: # remove cells with NaNs in covariates
                nan_idx = np.where(pd.isnull(cov_mat))
                nan_idx = np.unique(nan_idx[0])
                keep_idx = np.setdiff1d(list(range(dat.shape[0])),nan_idx)
                dat = dat[keep_idx,:]
                p_mat_pd = p_mat_pd[keep_idx,:]
                cov_mat = dat.obs[cov_names]
                log.info('Filtered {} cells with NaN in covariates'.format(len(nan_idx)))
        else:
            cov_mat = None

        # remove genes with 0 expression
        prev_count = dat.shape[1]
        scanpy.pp.filter_genes(dat, min_cells=1)
        after_count = dat.shape[1]
        log.info('Filtered {} genes with 0 expression across all cells'.format(prev_count - after_count))

        # remove cells with 0 counts
        prev_count = dat.shape[0]
        keep_cells = scanpy.pp.filter_cells(dat, min_counts=1, inplace=False)[0]
        dat = dat[keep_cells,:]
        p_mat_pd = p_mat_pd[keep_cells,:]
        if args.covariates is not None:
            cov_mat = cov_mat.loc[keep_cells,:]
        log.info('Filtered {} cells with 0 expression across all genes'.format(prev_count - np.sum(keep_cells)))

        log.info('Retained {} cells and {} genes for analysis'.format(dat.shape[0], dat.shape[1]))
        
        # normalize rows of p_mat if cell-pooled
        if overload_type == 'droplet':
            guides_per_cell = np.array(p_mat_pd.sum(axis = 1))
            vals = np.repeat(guides_per_cell, p_mat_pd.getnnz(axis=1))
            p_mat_pd.data = p_mat_pd.data / vals

        # convert raw counts to log(TP10K+1) 
        log.info('Converting raw counts to log(TP10K+1)...  ')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scanpy.pp.normalize_total(dat, target_sum = 10000)
        logmeanexp = np.squeeze(np.array(np.log(np.mean(dat.X, axis = 0)))) # for scaling later
        scanpy.pp.log1p(dat)

        # Running elastic net
        if args.elastic_net:
            log.info('Start running elastic net (--elastic-net set).')
            if args.elastic_net_genes_per_batch is not None:
                B = np.zeros((p_mat_pd.shape[1], dat.X.shape[1]))
                idxs = list(range(0, dat.X.shape[1], args.elastic_net_genes_per_batch)) + [dat.X.shape[1]]
                for i in range(len(idxs)-1):
                    log.info('Analyzing batch {} of {}.'.format(i+1, len(idxs)-1))
                    B_temp = ElasticNet(l1_ratio = 0.5, alpha = 0.0005, max_iter = 10000)
                    temp_mat = dat.X[:,idxs[i]:idxs[i+1]]
                    if scipy.sparse.issparse(temp_mat):
                        temp_mat = np.array(temp_mat.todense())
                    B_temp = B_temp.fit(p_mat_pd, temp_mat)
                    B_temp = B_temp.coef_.T
                    B[:,idxs[i]:idxs[i+1]] = B_temp
            else:
                B = ElasticNet(l1_ratio = 0.5, alpha = 0.0005, max_iter = 10000)
                temp_mat = dat.X
                if scipy.sparse.issparse(temp_mat):
                    temp_mat = np.array(temp_mat.todense())
                B = B.fit(p_mat_pd, temp_mat)
                B = B.coef_.T

        else:   
            # Running factorize recover
            if args.large_dataset:
                log.info('Start running FR-Perturb (--large-dataset set).')
                B, U, U_tilde, W, p_mat, alpha = factorize_recover_large_dataset(dat.X, p_mat_pd, pnames, args.rank, args.batches, control_perturbation_name=args.control_perturbation_name, cov_mat=cov_mat, log=log, alpha=args.lasso_alpha)        
            else:
                log.info('Start running FR-Perturb.')
    
                # Normalize matrix
                if cov_mat is not None:
                    log.info('Regressing out covariates...  ')
                    exp_mat = regress_covariates(dat.X, cov_mat)
                else:
                    exp_mat = dat.X
    
                if args.control_perturbation_name is not None:
                    log.info('Centering expression matrix based on control expression...  ')
                    # center expression matrix based on control expression
                    n_guides = p_mat_pd.sum(axis = 1)
                    nt_names = args.control_perturbation_name.split(',')
                    nt_idx = np.where(np.isin(pnames, nt_names))[0]
                    ctrl_idx = np.where(np.logical_and(n_guides == 1, p_mat_pd[:,nt_idx].sum(axis = 1) != 0))[0]
        
                    if len(ctrl_idx) < 100:
                        log.info('WARNING: Too few (<100) control cells, effect sizes will be computed relative to mean expression across all cells instead.')
                        ctrl_exp = exp_mat.mean(axis=0)
                    else:
                        ctrl_exp = exp_mat[ctrl_idx,:].mean(axis=0)
                else:
                    log.info('Centering expression matrix...  ')
                    ctrl_exp = exp_mat.mean(axis=0)
                
                exp_mat = exp_mat - ctrl_exp        
                exp_mat = np.asfortranarray(exp_mat.T)
                B, U, U_tilde, W, p_mat, alpha = factorize_recover(exp_mat, p_mat_pd, args.rank, args.spca_alpha, log=log, alpha=args.lasso_alpha)
    
            # Computing interactions
            if args.interaction_column_name is not None:
                log.info('Computing interaction terms...  ')
                resid = U_tilde - p_mat.dot(U)
    
                interaction_names = args.interaction_column_name.split(',')
                interaction_mat = dat.obs[interaction_names]
                interaction_mat = pd.get_dummies(interaction_mat)
                p_mat_dense = np.array(p_mat.todense())
    
                interaction_names = ['{}_{}'.format(pnames[x], interaction_mat.columns[y]) for x in range(p_mat_dense.shape[1]) for y in range(interaction_mat.shape[1])]
                
                pnames = np.concatenate([pnames, interaction_mat.columns.values, interaction_names])
                interaction_first_order_terms, interaction_second_order_terms = compute_interactions(p_mat_dense, interaction_mat.values, resid, alpha)
                U = np.r_[U, interaction_first_order_terms, interaction_second_order_terms]
                B = U.dot(W)
            
            log.info('Finished analysis.')

        # Scaling effects
        log.info('Scaling effects...  ')        
        B_scaled, _ = scale_effs(B, logmeanexp) 
        if args.gene_column_name is not None:
            gnames = dat.var[args.gene_column_name].values
        else:
            gnames = dat.var.index

        log.info('Outputting effects...  ')
        B_scaled = pd.DataFrame(data = np.transpose(B_scaled), index = gnames, columns = pnames)
        B_scaled = signif(B_scaled, 3)
        B_scaled.to_csv(args.out + '_LFCs.txt', sep = '\t')

        if args.temp_out:
            log.info('Outputting temporary files for batched permutation testing...  ')

            os.makedirs(os.path.dirname(args.temp_out), exist_ok=True)
            temp_out = [args, pnames, gnames, U_tilde, W, B, p_mat, alpha]
            if args.interaction_column_name:
                temp_out.append(interaction_mat)

            with open(args.temp_out + '.pkl', 'wb') as f:
                pickle.dump(temp_out, f)

        # Compute pvalues by permutation testing
        if args.compute_pval and not args.temp_out:
            if not args.fit_zero_pval:
                log.info('(--compute-pvals set) Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
                pvals = np.zeros((B.shape))
                for i in tqdm.tqdm(range(args.num_perms)):
                    perm_idxs = np.random.permutation(p_mat.shape[0])
                    p_mat_perm = p_mat[perm_idxs,:]
                    U_perm = Lasso(alpha = alpha, fit_intercept=False)
                    U_perm.fit(p_mat_perm, U_tilde)
                    U_perm = U_perm.coef_.T
                    B_perm = U_perm.dot(W)
                    if args.interaction_column_name: 
                        p_mat_perm = np.array(p_mat_perm.todense())
                        resid_perm = U_tilde - p_mat_perm.dot(U_perm)
                        interaction_mat_perm = interaction_mat.iloc[perm_idxs,:]
                        interaction_first_order_terms_perm, interaction_second_order_terms_perm = compute_interactions(p_mat_perm, interaction_mat_perm.values, resid_perm, alpha)                       
                        U_perm = np.r_[U_perm, interaction_first_order_terms_perm, interaction_second_order_terms_perm]
                        B_perm = U_perm.dot(W)
                    temp_indices = B < B_perm
                    pvals[temp_indices] = pvals[temp_indices] + 1
                pvals /= args.num_perms
                pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
                pvals *= 2 
                pvals = (pvals * args.num_perms + 1) / (args.num_perms + 1)
            else: 
                log.info('(--compute-pvals set) Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
                B_perms = np.empty((np.product(B.shape), args.num_perms))
                
                for i in tqdm.tqdm(range(args.num_perms)):
                    p_mat_perm = p_mat[np.random.permutation(p_mat.shape[0]),:]
                    U_perm = Lasso(alpha = alpha, fit_intercept=False)
                    U_perm.fit(p_mat_perm, U_tilde)
                    U_perm = U_perm.coef_.T
                    B_perm = U_perm.dot(W)
                    
                    if args.interaction_column_name: 
                        p_mat_perm = np.array(p_mat_perm.todense())
                        resid_perm = U_tilde - p_mat_perm.dot(U_perm)
                        interaction_mat_perm = interaction_mat.iloc[perm_idxs,:]
                        interaction_first_order_terms_perm, interaction_second_order_terms_perm = compute_interactions(p_mat_perm, interaction_mat_perm.values, resid_perm, alpha)                       
                        U_perm = np.r_[U_perm, interaction_first_order_terms_perm, interaction_second_order_terms_perm]
                        B_perm = U_perm.dot(W)
                    
                    B_perms[:,i] = np.ravel(B_perm)
                pvals = (B_perms < np.ravel(B)[:,np.newaxis]).sum(axis=1) / B_perms.shape[1]
                pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
                pvals *= 2 

                log.info('(--fit-zero-pval set) Fitting skew-normal distribution to effects with p=0 ({} total effects)...  '.format(np.sum(pvals == 0)))
                zero_indices = np.where(pvals == 0)[0]
                B_flattened = np.ravel(B)
                
                if args.multithreaded:
                    fitted_pvals = thread_map(lambda i: fit_skew_norm(B_flattened[i], B_perms[i,:]), zero_indices) 
                    
                else: 
                    fitted_pvals = []
                    for i in tqdm.tqdm(zero_indices):
                        t_star = B_flattened[i]
                        t_nulls = B_perms[i,:]
                        p = fit_skew_norm(t_star, t_nulls)
                        fitted_pvals.append(p)
                pvals[zero_indices] = fitted_pvals
                pvals = np.reshape(pvals, B.shape)
                log.info('Done.')

            qvals = sms.multitest.multipletests(pvals.flatten(), method = 'fdr_bh')[1]
            qvals = np.reshape(qvals, pvals.shape)
            pvals = pd.DataFrame(data = np.transpose(pvals), index = gnames, columns = pnames)
            qvals = pd.DataFrame(data = np.transpose(qvals), index = gnames, columns = pnames)
            pvals = signif(pvals, 3)
            qvals = signif(qvals, 3)
            pvals.to_csv(args.out + '_pvals.txt', sep = '\t')
            qvals.to_csv(args.out + '_qvals.txt', sep = '\t')

        # Running cross validation
        if args.cross_validate:
            log.info('(--cross-validate set) Running 2-fold cross validation... ')
            if args.large_dataset:
                args.cross_validate_runs = 1

            all_cor = []
            all_sc = []

            for i in range(args.cross_validate_runs):
                log.info('Run {} of {}... '.format(i+1, args.cross_validate_runs))
                n_cells = dat.X.shape[0]
                split1_idx = np.array(random.sample(range(n_cells), round(n_cells / 2)))
                split1_idx.sort()
                split2_idx = np.array(list(set(range(n_cells)) - set(split1_idx)))
                splits = [split1_idx, split2_idx]
                Bs = []
                for i, split_idx in enumerate(splits):
                    log.info('Split {}'.format(i+1))
                    if args.large_dataset:
                        B_temp,U_temp,U_tilde_temp,_,_,_ = factorize_recover_large_dataset(dat.X[split_idx,:], p_mat[split_idx,:], pnames, args.rank, args.batches, control_perturbation_name=args.control_perturbation_name, alpha=alpha, cov_mat=cov_mat.iloc[split_idx,:], log=log)
                    else:
                        B_temp,U_temp,U_tilde_temp,_,_,_ = factorize_recover(exp_mat[:,split_idx], p_mat[split_idx,:], args.rank, args.spca_alpha, log=log, alpha=alpha)

                    if args.interaction_column_name: 
                        p_mat_temp = np.array(p_mat[split_idx,:].todense())
                        resid_temp = U_tilde_temp - p_mat_temp.dot(U_temp)
                        interaction_mat_temp = interaction_mat.iloc[split_idx,:]
                        interaction_first_order_terms_temp, interaction_second_order_terms_temp = compute_interactions(p_mat_temp, interaction_mat_temp.values, resid_temp, alpha)                       
                        U_temp = np.r_[U_temp, interaction_first_order_terms_temp, interaction_second_order_terms_temp]
                        B_temp = U_temp.dot(W)
                    B_temp,_ = scale_effs(B_temp, logmeanexp)
                    Bs.append(B_temp)

                temp_cors = []
                temp_sc = []

                # separately cross-validate for interaction effects
                if args.interaction_column_name:
                    if args.compute_pval and not args.temp_out:
                        cnames = ['{}_{}'.format(x, y) for x in ['First_order_perturbation', 'First_order_condition', 'Second_order_interaction'] for y in ['top_0.1%', 'top_1%', 'q<0.05', 'q<0.2']]
                    else:
                        cnames = ['{}_top_{}%'.format(x, y) for x in ['First_order_perturbation', 'First_order_condition', 'Second_order_interaction'] for y in [0.1, 1]]
                    etype_idxs = [0] + list(np.cumsum([p_mat.shape[1], interaction_mat.shape[1]])) + [B.shape[0] + 1]
                    for i in range(3):
                        temp_Bs = [x[etype_idxs[i]:etype_idxs[i+1],:].flatten() for x in Bs]
                        sort_idxs = np.argsort(-np.abs(B_scaled.iloc[:,etype_idxs[i]:etype_idxs[i+1]].values.flatten()))
                        cutoffs = np.round(np.array([0.001, 0.01]) * len(sort_idxs)).astype(int)

                        for cutoff in cutoffs:
                            cor = np.corrcoef(temp_Bs[0][sort_idxs[:cutoff]], temp_Bs[1][sort_idxs[:cutoff]])[0,1]
                            sc = sign_concord(temp_Bs[0][sort_idxs[:cutoff]], temp_Bs[1][sort_idxs[:cutoff]])
                            temp_cors.append(cor)
                            temp_sc.append(sc)

                        if args.compute_pval and not args.temp_out: # also compute correlation using significant effects
                            temp_qvals = qvals.iloc[:,etype_idxs[i]:etype_idxs[i+1]]

                            for cutoff in [0.05, 0.2]:
                                idxs = (temp_qvals.values < cutoff).flatten()
                                if np.sum(idxs) == 0:
                                    cor = np.nan
                                    sc = np.nan
                                else:    
                                    cor = np.corrcoef(temp_Bs[0][idxs], temp_Bs[1][idxs])[0,1]
                                    sc = sign_concord(temp_Bs[0][idxs], temp_Bs[1][idxs])
                                temp_cors.append(cor)
                                temp_sc.append(sc)
                        
                else:
                    cnames = ['Top_{}%_effects'.format(x) for x in [0.1, 1]]
                    Bs = [x.flatten() for x in Bs]
                    sort_idxs = np.argsort(-np.abs(B_scaled.values.flatten()))
                    cutoffs = np.round(np.array([0.001, 0.01]) * len(sort_idxs)).astype(int)
    
                    for cutoff in cutoffs:
                        cor = np.corrcoef(Bs[0][sort_idxs[:cutoff]], Bs[1][sort_idxs[:cutoff]])[0,1]
                        sc = sign_concord(Bs[0][sort_idxs[:cutoff]], Bs[1][sort_idxs[:cutoff]])
                        temp_cors.append(cor)
                        temp_sc.append(sc)
    
                    if args.compute_pval and not args.temp_out: # also compute correlation using significant effects
                        cnames.extend(['q<{}_effects'.format(x) for x in [0.05, 0.2]])
                        
                        for cutoff in [0.05, 0.2]:
                            idxs = (qvals.values < cutoff).flatten()
                            if np.sum(idxs) == 0:
                                cor = np.nan
                                sc = np.nan
                            else:    
                                cor = np.corrcoef(Bs[0][idxs], Bs[1][idxs])[0,1]
                                sc = sign_concord(Bs[0][idxs], Bs[1][idxs])
                            temp_cors.append(cor)
                            temp_sc.append(sc)

                all_cor.append(temp_cors)
                all_sc.append(temp_sc)
                    
            all_cor = np.array(all_cor).mean(axis=0)
            all_sc = np.array(all_sc).mean(axis=0)
            cv_out = pd.DataFrame([all_cor, all_sc], columns=cnames, index=['Correlation', 'Sign_concordance']).T
            cv_out.to_csv(args.out + '_twofold_cross_validation_results.txt', sep = '\t')

        log.info('All done!')       
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.info(traceback.format_exc(ex))
        raise
    finally:
        log.info('Analysis finished at {T}'.format(T=time.ctime()))
        time_elapsed = round(time.time() - start_time, 2)
        log.info('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

