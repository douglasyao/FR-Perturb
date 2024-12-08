#!/usr/bin/env python
# coding: utf-8

# In[2]:


import time, sys, traceback, argparse
import os
import random
import warnings
import pickle
import tqdm
import glob
from sklearn.linear_model import Lasso, LassoCV
from FR_Perturb.FR_Perturb import *
from FR_Perturb.utils import *


# In[4]:


parser = argparse.ArgumentParser()

### Args
parser.add_argument('--input-prefix', default=None, type=str,
                    help='Prefix that should match --temp-out flag in main script.')
parser.add_argument('--batch-number', default=None, type=int,
                    help='Batch number')
parser.add_argument('--perms-per-batch', default=25, type=int,
                    help='How many permutations to perform per batch.')
parser.add_argument('--combine-pval-files', default=False, action='store_true',
                    help='After permutations are computed, set this flag to combine them to compute p-values')
parser.add_argument('--out', default=None, type=str,
                    help='Output file prefix')
                    
# input_args = parser.parse_args(['--input-prefix', '/n/scratch/users/d/dwy6/scperturb/temp_files/asdf',
#                   '--batch-number', '1'])
input_args = parser.parse_args()

if not input_args.input_prefix:
    raise ValueError('Must specify --input-prefix')

if input_args.combine_pval_files and not input_args.out:
    raise ValueError('If --combine-pval-files is set, must specify --out')
    
with open(input_args.input_prefix + '.pkl', 'rb') as f:
    dat = pickle.load(f)

args = dat[0]
pnames = dat[1]
gnames = dat[2]
U_tilde = dat[3]
W = dat[4]
B = dat[5]
p_mat = dat[6]
alpha = dat[7]
if args.interaction_column_name:
    interaction_mat = dat[8]

if input_args.combine_pval_files:
    all_files = glob.glob(input_args.input_prefix + '_batch_*_perms.pkl')
    print('Combining permutation files...')
    print('Permutations per batch: {}'.format(input_args.perms_per_batch))
    print('Loading {} batch files'.format(len(all_files)))

    total_perms = len(all_files) * input_args.perms_per_batch
    print('Total permutations: {}'.format(total_perms))

    pvals = np.zeros((B.shape))
    for file in all_files:
        with open(file, 'rb') as f:
            temp = pickle.load(f)
        pvals = pvals + temp
    
    pvals /= total_perms
    pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
    pvals *= 2 
    pvals = (pvals * total_perms + 1) / (total_perms + 1)

    qvals = sms.multitest.multipletests(pvals.flatten(), method = 'fdr_bh')[1]
    qvals = np.reshape(qvals, pvals.shape)
    pvals = pd.DataFrame(data = np.transpose(pvals), index = gnames, columns = pnames)
    qvals = pd.DataFrame(data = np.transpose(qvals), index = gnames, columns = pnames)
    pvals = signif(pvals, 3)
    qvals = signif(qvals, 3)
    pvals.to_csv(input_args.out + '_pvals.txt', sep = '\t')
    qvals.to_csv(input_args.out + '_qvals.txt', sep = '\t')
    print('Done!')
    
    sys.exit()
    
np.random.seed(input_args.batch_number)
pvals = np.zeros((B.shape))
print('Running permutations... ')
for i in tqdm.tqdm(range(input_args.perms_per_batch)):
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

print('Done!')
with open(input_args.input_prefix + '_batch_{}_perms.pkl'.format(input_args.batch_number), 'wb') as f:
    pickle.dump(pvals, f)



