#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy
import time, sys, traceback, argparse
import os
import tqdm
import random
import logging
from tqdm.contrib.concurrent import thread_map
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.linear_model import Lasso
from FR_Perturb.FR_Perturb import *
from FR_Perturb.utils import *


# In[3]:


os.environ['KMP_WARNINGS'] = 'off'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* Factorize-Recover for Perturb-seq analysis (FR-Perturb)\n"
MASTHEAD += "* by Douglas Yao 2022 \n"
MASTHEAD += "*********************************************************************\n"

parser = argparse.ArgumentParser()

### Required flags
parser.add_argument('--input-h5ad', default=None, type=str,
                    help='h5ad file (from the AnnData package) containing raw gene expression counts for all cells')
parser.add_argument('--control-perturbation-name', default=None, type=str,
                    help='Comma-separated list of perturbation names that represent control perturbations')
parser.add_argument('--out', default=None, type=str,
                    help="Output prefix (including directory) for effect sizes")

### Either need to specify 
parser.add_argument('--input-perturbation-matrix', default=None, type=str,
                    help='Whitespace-delimited file containing a table with columns corresponding to cells and rows corresponding to perturbations. Cells containing a given perturbation should be indicated with a "1", otherwise "0".')
### Or
parser.add_argument('--perturbation-column-name', default=None, type=str,
                    help='Column name in the cell annotation matrix (stored in `.obs` in the h5ad object) that corresponds to which perturbations are in the cell.')
parser.add_argument('--perturbation-delimiter', default='--', type=str,
                    help='Delimiter used to separate perturbations in each cell (if the cell contains multiple perturbations). If the delimiter includes a double-dash --, then please enclose in both double and single quotes (e.g. '"--"'). This is a hack, sorry.')

### Optional
parser.add_argument('--separate-control-perturbations', default=False, action='store_true',
                    help='Separately infer effects for each individual control perturbation. By default, all separate control perturbations are combined into a single "control" category.')
parser.add_argument('--compute-pval', default=False, action='store_true',
                    help='Whether or not to compute p-values for all effect size estimates by permutation testing (WARNING: this can take a long time).')           
parser.add_argument('--rank', default=20, type=int,
                    help='Hyperparameter determining the rank of the matrix during the factorize step')
parser.add_argument('--lambda1', default=0.1, type=float,
                    help='Hyperparameter determining the sparsity of the factor matrix during the factorize step of the method. Higher value = more sparse.')
parser.add_argument('--lambda2', default=10, type=float,
                    help='Hyperparameter determining the sparsity of learned effects during the recover step of the method. Higher value = more sparse.')
parser.add_argument('--covariates', default=None, type=str,
                    help='Comma-separated list of covariate names to regress out of the expression matrix (names must match the column names in the meta-data of the h5ad object)')
parser.add_argument('--large-dataset', default=False, action='store_true',
                    help='Set this flag if dataset is too large to load into memory (>1M cells). Uses memory-efficient randomized partial SVD during factorize step of factorize-recover. --lambda1 parameter will be ignored.')
parser.add_argument('--batches', default=5, type=int,
                    help='Batch size when --large-dataset is set. Increasing this number reduces the amount of memory but increases running time.')
parser.add_argument('--guide-pooled', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from guide pooling')
parser.add_argument('--cell-pooled', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from cell pooling')
parser.add_argument('--num-perms', default=1000, type=int,
                    help='Number of permutations when doing permutation testing')
parser.add_argument('--fit-zero-pval', default=False, action='store_true',
                    help='Compute p-values by fitting skew-normal distribution to null distribution (allows for p-values below 1/num_perms, but significantly increases compute time)')
parser.add_argument('--multithreaded', default=False, action='store_true',
                    help='Use multithreading to fit skew-normal distributions, which can substantially reduce compute time')
parser.add_argument('--output-factor-matrices', default=False, action='store_true',
                    help='Whether or not to output the latent gene expression factor matrices in addition to the full effect sizes matrix')


## Cross-validation                    
parser.add_argument('--cross-validate', default=False, action='store_true',
                    help='Perform 2-fold cross-validation on input data')
parser.add_argument('--cross-validate-runs', default=2, type=int,
                    help='Number of times to repeat cross-validation. Cross-validation accuracy is averaged across runs.')


# In[ ]:


if __name__ == '__main__':

    args = parser.parse_args()

    # args = parser.parse_args(['--input-h5ad', '/n/scratch/users/d/dwy6/frperturb_testing/GSM6858448_KO_cell_pooled.h5ad',
    #                   '--perturbation-column-name', 'Guides_collapsed_by_gene',
    #                   '--perturbation-delimiter', '"--"',
    #                   '--control-perturbation-name', 'non-targeting,safe-targeting',
    #                   '--out', '/n/scratch/users/d/dwy6/asdf',
    #                   '--covariates', 'Total_RNA_count,Percent_mitochondrial_reads,S_score,G2M_score'])

    # args = parser.parse_args(['--input-h5ad', '/n/scratch/users/d/dwy6/scperturb/data/GasperiniShendure2019_atscale.h5ad',
    #                   '--perturbation-column-name', 'barcode',
    #                   '--control-perturbation-name', 'AATGAGGAGCAAACGAAAAT,ACGAAATGTTTCATGACCAA,ATAGATTTACGTTACTCTCT,ATTAGCATCAGGTAGACTAA,CCATAAAGAATTCGGTGTAG,CGAAAATGGGTGCGTGGACT,CGAAGGATCAAAGCGACTTT,CGTATTTCTCTTAGAGATCA,GAGAAAAGTCGTTCTACAAA,GCGTCCCTAGGGTATTTTAC,TCTTGGCCTGCTTGGTGTCT,CTCACCTGTGGGAGTAACGG,ACCGCTGCCCTAACCACTGG,ACTGCCTCGCGATTGACTGG,AGAGTGTAGATACGTTGGTG,AGCCTAACGATCGGACCGAG,AGCTCTTGCGGTAAATCGTG,ATGCGGCGAACTGCCGTAAG,CAAAGGCTGAGGAGCGGTGG,CAATCAGCGACTGCTGCTGG,CACCACTGGCATGCTTGAAG,CACCGGGACTCGGGAGTGCG,CAGGATCGCTATCAGCACGG,CCACATACGAAATGCCAACG,CCACCTTCGAAGTCCGTATG,CCCTCAAGGGCACGTTTAGG,CGGTGGGGCCGCTCACTAGG,CTATTCGCCAACTTGTTACG,CTCACCGTAGAGACGGCACG,CTGAACTTCGGTAAAGGCGG,CTGCAACGTGGTGCTGCCGG,GAACTCAATAGTGAGAGCGG,GACCTCCTGTGATCAGGTGG,GCGAACAGCGAGGGCCCCTG,GCGACATTTGGGTCGCGAAG,GCGGCGCCTGTAAGCCCCTG,GCTGTATATCGGCGCCCCGG,GGACGAGTAACCTGCCGGGG,GGCATCGGGCAGACGGGATG,GGGCAGTCGTAACTCTAAGG,GGGGACTTGCTCCAGGACGG,GTACCCGCTAAGGAGCGCGG,GTGCCCGGAGCTACATGGCG,TAAACTCATGGTCCTGAAGG,TAAGGGCGACGTTCGGTAAG,TAAGGGCGCCCGAATCGCCG,TAGCTTGCACTCCCCTTGTG,TATTCGTACCGGGCAGCAGG,TGTCCAGTTCCGTAGGATGG,TGTTGTGAGAGCATCCGGAG,TTGGGAATGCGAGGCAAAGG,CTAAAGCATTGGCTGAGAAG,AATCATGGTGGAAGGTGAAG,GATATGTTATTGAAATCTAA,TCTAATCTCAGCTACTTGGG,GTAGTTCACATAATCCCTGT,GCCTTGTTTCTATTTGTAGA,AACACAACACACCAAAACTG,GATGTTTTATATAATCCCCA,AAGTTGACTCTACATAGCAG,AATATTCTCCCTCATTCTGG,AATCCTCTAATGGACGAAGA,ACACACACAGAGAAATATAG,AGAATTTCTTTGAACCCGGG,AGAGGTAACCAAAATAGCAA,AGATACCTATGGCCATATAG,AGGAATCCAGCAATCCCACA,AGTTCTAGAGCACTGAGCAA,ATAGAAGGAGGATGAGGGGA,ATAGGCCAGTTACAATTTGG,ATATGTAACCTCCAGAATGA,ATATTCAGCAGCTAAAGCAT,ATGTGTATATTAACATAGGG,ATTGGCAGTCTCTAAGAAGT,ATTGGTATCCGTATAAGCAG,CACCCTCCTGGGGAAGACCA,CCTGACAATCAATCTCCAAG,CGAGTTGTAAGCCCTTAAAA,CTATAGGGTATGAAGAGCAA,GATGAAGAAGTATAAAGCAG,GCAAATGCTTCATCACCCCA,GCCTGCAGAAAGCTTCCCTG,GCTCTAATGAACAGAATGGG,GCTGGGGAGACCCAACCCAG,GGGCAACCCCAGGAAGACCG,GTAGCCTCTGTTCCTCAGTA,GTAGCTTTATGAGAACCACT,GTGGCCCTACTGAGGAACAG,TACAACTGCATTACATGCCA,TATCTCAGAGTACTATTCCA,TCAGGGGTCGATCTTTAACC,TCCGCAGTCAAAAGACCGAG,TCTGCTAAACTGCCTACACA,TGAACAATACTCCAGTACAT,TGCATCTTCTGAAATGGCAG,TTAAAATTGATTCTGCCACT,TTAAGGGCTTACAACTCGAA,TTAATTCCTCTGGCGCCGCT,TTAGCTTTGAAAACACACAC,TTATCTCTATTTGACAGACG,TTTCTCTGTAAGTTACCATG',
    #                   '--covariates', 'total_umis,nperts,percent.mito,prep_batch',
    #                   '--perturbation-delimiter', '_',
    #                   '--out', '/n/scratch/users/d/dwy6/asdf'])

    logging.basicConfig(
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.FileHandler(args.out + ".log"),logging.StreamHandler()],
        level=logging.INFO,
    )
    log = logging.getLogger(__name__)

    try:
        if args.large_dataset:
            args.lambda2 = 0.0001

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
        if not args.control_perturbation_name:
            raise ValueError('Must specify --control-perturbation-name')
        if not args.out:
            raise ValueError('Must specify --out')
        if args.guide_pooled and args.cell_pooled:
            raise ValueError('Only one of --guide-pooled and --cell-pooled should be set')
        elif args.cell_pooled:
            overload_type = 'droplet'
        else:
            overload_type = 'guide' # use this version by default

        log.info('Loading input data from {}...  '.format(args.input_h5ad))
        dat = scanpy.read_h5ad(args.input_h5ad)
        
        # get perturbation matrix
        if args.input_perturbation_matrix:
            p_mat_pd = pd.read_csv(args.input_perturbation_matrix, index_col = 0, delim_whitespace=True)
            if not dat.obs.index.equals(p_mat_pd.columns):
                raise ValueError('Cell names in perturbation matrix do not match cell names in expression matrix')
        else:
            perts = dat.obs[args.perturbation_column_name]
            if args.perturbation_delimiter[0] == '\"' and args.perturbation_delimiter[-1] == '\"': # if delimiter includes double dash
                args.perturbation_delimiter = args.perturbation_delimiter[1:-1]

            mlb = MultiLabelBinarizer(sparse_output=True)
            p_mat_pd = mlb.fit_transform(perts.str.split(args.perturbation_delimiter))
            pnames=mlb.classes_

        pert_idx = np.where(np.isin(pnames, args.control_perturbation_name.split(',')))[0]
        if len(pert_idx) != len(args.control_perturbation_name.split(',')):
            raise ValueError('Provided control perturbation names not found in data. Please also double check that delimiter is specified correctly.')

        if not args.separate_control_perturbations and len(args.control_perturbation_name.split(',')) > 1: # combine control perturbations into one category
            ctrl_perts = p_mat_pd[:,pert_idx]
            ctrl_perts = ctrl_perts.max(axis = 1)
            mask = np.ones(p_mat_pd.shape[1], dtype=bool)
            mask[pert_idx] = False
            p_mat_pd = scipy.sparse.hstack([p_mat_pd[:,mask], ctrl_perts], format='csr')
            pnames = np.append(pnames[mask], 'control')
            args.control_perturbation_name = 'control'
        
        # center rows of p_mat
        if overload_type == 'droplet':
            guides_per_cell = np.array(p_mat_pd.sum(axis = 1))
            vals = np.repeat(guides_per_cell, p_mat_pd.getnnz(axis=1))
            p_mat_pd.data = p_mat_pd.data / vals

        # get covariates
        if args.covariates is not None:
            cov_names = args.covariates.split(',')
            cov_mat = dat.obs[cov_names]
        else:
            cov_mat = None

        # convert raw counts to log(TP10K+1) 
        log.info('Converting raw counts to log(TP10K+1)...  ')
        scanpy.pp.normalize_total(dat, target_sum = 10000)
        logmeanexp = np.squeeze(np.array(np.log(np.mean(dat.X, axis = 0)))) # for scaling later
        scanpy.pp.log1p(dat)
        
        # Running factorize recover
        if args.large_dataset:
            log.info('Start running FR-Perturb (--large_dataset set).')
            B, U, U_tilde, W, p_mat = factorize_recover_large_dataset(dat.X, p_mat_pd, pnames, args.rank, args.lambda2, args.batches, args.control_perturbation_name, cov_mat=cov_mat, log=log)
        else:
            log.info('Start running FR-Perturb.')
            p_mat = pd.DataFrame.sparse.from_spmatrix(p_mat_pd, columns=pnames)

            # Normalize matrix
            if cov_mat is not None:
                log.info('Regressing out covariates...  ')
                exp_mat = regress_covariates(dat.X, cov_mat)
        
            log.info('Centering expression matrix based on control expression...  ')
            # center expression matrix based on control expression
            n_guides = p_mat.values.sum(axis = 1)
            nt_names = args.control_perturbation_name.split(',')
            ctrl_idx = np.logical_and(n_guides == 1, p_mat.loc[:,nt_names].sum(axis = 1).values != 0)
            
            if np.sum(ctrl_idx) < 100:
                log.info('WARNING: Too few (<100) control cells, effect sizes will be computed relative to mean expression across all cells instead.')
            else:
                ctrl_exp = exp_mat[ctrl_idx,:].mean(axis = 0)
                exp_mat = exp_mat - ctrl_exp
            exp_mat = np.asfortranarray(exp_mat.T)
            B, U, U_tilde, W, p_mat = factorize_recover(exp_mat, p_mat, args.rank, args.lambda1, args.lambda2, log=log)

        log.info('Finished running FR-Perturb.')

        # Scaling effects
        log.info('Scaling effects...  ')
        B_scaled, _ = scale_effs(B, logmeanexp)            

        # Compute pvalues by permutation testing
        if args.compute_pval:
            if not args.fit_zero_pval:
                log.info('Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
                pvals = np.zeros((B.shape))
                for i in tqdm.tqdm(range(args.num_perms)):
                    p_mat_perm = p_mat[np.random.permutation(p_mat.shape[0]),:]
                    if args.large_dataset:
                        U_perm = Lasso(alpha = args.lambda2, fit_intercept=False)
                        U_perm.fit(p_mat_perm, U_tilde)
                        U_perm = U_perm.coef_.T
                    else:
                        U_perm = spams.lasso(U_tilde, D=np.asfortranarray(p_mat_perm), lambda1=args.lambda2, verbose=False)
                    B_perm = U_perm.dot(W)
                    temp_indices = B < B_perm
                    pvals[temp_indices] = pvals[temp_indices] + 1
                pvals /= args.num_perms
                pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
                pvals *= 2 
                pvals = (pvals * args.num_perms + 1) / (args.num_perms + 1)
            else:
                log.info('Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
                B_perms = np.empty((np.product(B.shape), args.num_perms))
                
                for i in tqdm.tqdm(range(args.num_perms)):
                    p_mat_perm = p_mat[np.random.permutation(p_mat.shape[0]),:]
                    if args.large_dataset:
                        U_perm = Lasso(alpha = args.lambda2, fit_intercept=False)
                        U_perm.fit(p_mat_perm, U_tilde)
                        U_perm = U_perm.coef_.T
                    else:
                        U_perm = spams.lasso(U_tilde, D=np.asfortranarray(p_mat_perm), lambda1=args.lambda2, verbose=False)
                    B_perm = U_perm.dot(W)
                    B_perms[:,i] = np.ravel(B_perm)
                pvals = (B_perms < np.ravel(B)[:,np.newaxis]).sum(axis=1) / B_perms.shape[1]
                pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
                pvals *= 2 

                log.info('Fitting skew-normal distribution to effects with p=0 ({} total effects)...  '.format(np.sum(pvals == 0)))
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
            pvals = pd.DataFrame(data = np.transpose(pvals), index = dat.var.index, columns = pnames)
            qvals = pd.DataFrame(data = np.transpose(qvals), index = dat.var.index, columns = pnames)
            pvals = signif(pvals, 3)
            qvals = signif(qvals, 3)

        log.info('Outputting results...  ')
        B_scaled = pd.DataFrame(data = np.transpose(B_scaled), index = dat.var.index, columns = pnames)
        B_scaled = signif(B_scaled, 3)
        B_scaled.to_csv(args.out + '_LFCs.txt', sep = '\t')
        if args.compute_pval:
            pvals.to_csv(args.out + '_pvals.txt', sep = '\t')
            qvals.to_csv(args.out + '_qvals.txt', sep = '\t')

        # Running cross validation
        if args.cross_validate:
            log.info('Running 2-fold cross validation... ')
            if args.large_dataset:
                args.cross_validate_runs = 1
            cutoffs=[0.001, 0.01]
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
                        B_temp,_,_,_,_ = factorize_recover_large_dataset(dat.X[split_idx,:], p_mat[split_idx,:], pnames, args.rank, args.lambda2, args.batches, args.control_perturbation_name, cov_mat=cov_mat.iloc[split_idx,:], log=log)
                    else:
                        B_temp,_,_,_,_ = factorize_recover(exp_mat[:,split_idx], p_mat[split_idx,:], args.rank, args.lambda1, args.lambda2, log=log)
                    B_temp,_ = scale_effs(B_temp, logmeanexp)
                    Bs.append(B_temp)
                cor, sc = compute_correlation(Bs[0], Bs[1], B_scaled.values, cutoffs)
                all_cor.append(cor)
                all_sc.append(sc)
            all_cor = np.array(all_cor).mean(axis=0)
            all_sc = np.array(all_sc).mean(axis=0)
            cv_out = pd.DataFrame([all_cor, all_sc], columns=['Top_{}%_effects'.format(x*100) for x in cutoffs], index=['Correlation', 'Sign_concordance']).T
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


