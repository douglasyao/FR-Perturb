# Factorize-Recover for Perturb-seq analysis (FR-Perturb)
FR-Perturb is a command-line tool for estimating effect sizes on gene expression from Perturb-seq data, based on the factorize-recover algorithm from [Sharan et al. 2019 ICML](http://proceedings.mlr.press/v97/sharan19a/sharan19a-supp.pdf). The code is still in preliminary form and will be updated in the future. The method is written in both R (used to run the analysis in the paper) and Python. We recommend using the Python version because it is much faster, though the output may be slightly different than the R version due to the use of different packages to perform inference. 

To use, clone the repository:
```
git clone https://github.com/douglasyao/FR-Perturb.git
cd FR-Perturb
```

Next, activate the conda environment:

```
conda env create --file environment.yml
conda activate FR-Perturb
```

FR-Perturb is now ready to use. For comments/questions/issues with running the code, please contact douglasyao@g.harvard.edu. 

# Usage

The following command should be used to estimate effect sizes:

`
./run_FR_Perturb.py --input-h5ad [INPUT_H5AD] --input-perturbation-matrix [INPUT_PERTURBATION_MATRIX] --control-perturbation-name [CONTROL_PERTURBATION_NAME] --covariates [COVARIATES] --compute-pval --fit-zero-pval --multithreaded --out [OUT]
`

### Required input
**1**. `--input-h5ad [INPUT_H5AD]` specifies an h5ad file (from the AnnData package) containing raw gene expression counts for all cells. Sample-level covariates should be found in the `obs` variable (see [here](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) for info on AnnData objects).   

**2**. `--input-perturbation-matrix [INPUT_PERTURBATION_MATRIX]` specifies a whitespace-delimited file containing a table with columns corresponding to cells and rows corresponding to perturbations. Cells containing a given perturbation should be indicated with a "1", otherwise "0". The column names should match the cell names in the h5ad file i.e. the row names of the `obs` variable (otherwise the method will throw an error). The row names should correspond to the perturbation names. An example is shown below:

```
CELL_1  CELL_2 CELL_3
PERT_1  1 0 0 
PERT_2  0 1 0
PERT_3  0 0 1
```

**3**. `--control-perturbation-name [CONTROL_PERTURBATION_NAME]` specifies the perturbation name corresponding to a control guide. The name must be present in the row names of the input perturbation matrix. If multiple rows correspond to control guides, then their names can be provided as a list separated by commas, e.g. `ctrl_pert1,ctrl_pert2,ctrl_pert3`.

**4**. `--out [OUT]` specifies the output file prefix, including the directory. The output name will be the prefix with `_LFCs.txt` appended to it. 

### Optional input

**5**. `--covariates [COVARIATES]` specifies a comma-separated list of covariate names to regress out of the expression matrix, e.g. `cov_1,cov_2,cov_3`. These names must be present in the column names of the `obs` variable in the h5ad file. 

**6**. `--compute-pval` specifies whether or not to compute p-values for all effect sizes by permutation testing.

**7**. `--fit-zero-pval` specifies whether or not to fit a skew-normal distribution the null distribution of effect sizes with p=0 from permutation testing. We recommend setting this flag as it allows for p-values below 1/num_perms, though it substantially increases compute time. 

**8**. `--multithreaded` specifies whether or not use multithreading when fitting skew-normal distributions (can substantially reduce compute time). 

**9**. `--num-perms [NUM_PERMS]` specifies the number of permutations to perform for permutation testing (default 10,000; 500 if `--fit-zero-pval` is set). 

**10**. `--rank [RANK]` specifies a hyperparameter determining the rank of the matrix during the factorize step (default 20).

**11**. `--lambda1 [LAMBDA1]` specifies a hyperparameter determining the sparsity of the left factor matrix during the factorize step; higher numbers = more sparse. Default 0.1.

**12**. `--lambda2 [LAMBDA2]` specifies a hyperparameter determining the sparsity of the effect sizes on latent factors during the recovery step; higher numbers = more sparse. Default 10. 

**13**. `--guide-pooled` specifies whether to run the version of FR-Perturb that assumes data is generated from guide pooling. 

**14**. `--cell-pooled` specifies whether to run the version of FR-Perturb that assumes data is generated from cell pooling. 

# Output

FR-Perturb will output a whitespace-delimited text file containing the log-fold changes of expression relative to expression in cells containing control guides, with rows corresponding to genes and columns corresponding to perturbations. If `--compute-pval` is set, the FR-Perturb will also output a file containing p-values (with the same rows and columns as the effect size file) and a file containing q-values (computed using the Benjamini-Hochberg procedure). An example is shown below:

```
PERT_1  PERT_2  PERT_3
GENE_1  0.39 -0.9 2.1 
GENE_2  0.32 -0.82 3.1
GENE_3  -0.61 0.24 1.2
```



