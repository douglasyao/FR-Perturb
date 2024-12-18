# Factorize-Recover for Perturb-seq analysis (FR-Perturb)
FR-Perturb is a command-line tool for estimating effect sizes on gene expression from Perturb-seq data, based on the factorize-recover algorithm from [Sharan et al. 2019 ICML](http://proceedings.mlr.press/v97/sharan19a/sharan19a-supp.pdf). See our publication ["Scalable genetic screening for regulatory circuits using compressed Perturb-seq" (Yao et al. 2023 Nature Biotechnology)](https://www.nature.com/articles/s41587-023-01964-9) for more details. The method is written in both R (used to run the analysis in the paper) and Python. We recommend using the Python version because it is much faster, though the output may be slightly different than the R version due to the use of different packages to perform inference. 

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

FR-Perturb is now ready to use. For comments/questions/issues with running the code, please contact douglasyao2@gmail.com. 

# Table of contents
- [Input](#input)
- [Output](#output)
- [Options for large datasets](#options-for-large-datasets)
- [Other options](#other-options)
  - [Gene names](#gene-names)
  - [Significance testing parameters](#significance-testing-parameters)
  - [Model hyperparameters](#model-hyperparameters)
  - [Running other perturbation effect inference methods](#running-other-perturbation-effect-inference-methods)
  - [Cross-validation](#cross-validation)

# Input

In order to estimate perturbation effect sizes, two things need to provided:

**1**. A raw expression count matrix, stored as an .h5ad file (from the AnnData package).

**2**. The perturbation status of each cell. This data can be specified in two ways:  
- **Option 1**. A cell x perturbation binary indicator matrix. This is provided as a separate whitespace-delimited file whose contents should look something like this:
```
Cell_Barcode  Pert_1  Pert_2 Control  ...
Cell_1  1 0 0 
Cell_2  0 0 1 
Cell_3  0 1 0 
Cell_4  1 1 0 
...
```
- **Option 2**. Info stored in the meta-data of the expression matrix (i.e. the `obs` variable of the AnnData object). When the expression matrix is loaded into Python as `data`, then `data.obs` should contain look something like this:
```
Cell_Barcode  Perturbation  Covariate_1  Covariate_2  ...
Cell_1  Pert_1  0.04  45
Cell_2  Control  0.14  69
Cell_3  Pert_2  0.29  12
Cell_4  Pert_1:Pert_2  0.22  45
...
```
### Option 1

For Option 1, the following command should be used to estimate effect sizes:

`
./run_FR_Perturb.py --input-h5ad [INPUT_H5AD] --input-perturbation-matrix [INPUT_PERTURBATION_MATRIX] --control-perturbation-name [CONTROL_PERTURBATION_NAME] --covariates [COVARIATES] --compute-pval --out [OUT]
`

**1**. `--input-h5ad [INPUT_H5AD]` specifies an h5ad file name (from the AnnData package) containing raw gene expression counts for all cells. Sample-level covariates should be found in the `obs` variable (see [here](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) for info on AnnData objects).   

**2**. `--input-perturbation-matrix [INPUT_PERTURBATION_MATRIX]` specifies a file name containing a whitespace-delimited table with columns corresponding to perturbations and rows corresponding to cells. Cells containing a given perturbation should be indicated with a "1", otherwise "0". The row names in this file should match the cell names in the h5ad file, i.e. the row names of the `obs` variable (otherwise the method will throw an error). The column names should correspond to the perturbation names. An example is shown below:

```
Cell_Barcode  Pert_1  Pert_2 Control  ...
Cell_1  1 0 0 
Cell_2  0 0 1 
Cell_3  0 1 0 
Cell_4  1 1 0 
...
```

**3**. `--control-perturbation-name [CONTROL_PERTURBATION_NAME]` (Optional). Typically, perturbation effect sizes are computed relative to cells containing a control perturbation, so these control cells should be present in the expression matrix. This flag specifies the perturbation name corresponding to a control guide. The name must be present in the column names of the input perturbation matrix. If multiple different columns correspond to control guides, then their names can be provided as a list separated by commas, e.g. `ctrl_pert1,ctrl_pert2,ctrl_pert3`. If this flag is not set, then perturbation effect sizes will be computed relative to the mean expression across all cells. 

**4**. `--covariates [COVARIATES]` (Optional). This optional flag specifies a comma-separated list of covariate names to regress out of the expression matrix, e.g. `cov_1,cov_2,cov_3`. These names must be present in the column names of the `obs` variable in the h5ad file. 

**5**. `--compute-pval` (Optional). This optional flag specifies whether or not to compute p-values for all effect sizes by permutation testing. If set, two additional files with suffixes `_pvals.txt` and `_qvals.txt` will be output. 

**6**. `--out [OUT]` specifies the output file prefix, including the directory. The output name for the effect sizes will be the prefix with `_LFCs.txt` appended to it. 
<br /><br />
### Option 2

For Option 2, the following command should be used to estimate effect sizes:

`
./run_FR_Perturb.py --input-h5ad [INPUT_H5AD] --perturbation-column-name [INPUT_PERTURBATION_MATRIX] --perturbation-delimiter [PERTURATION_DELIMITER] --control-perturbation-name [CONTROL_PERTURBATION_NAME] --covariates [COVARIATES] --compute-pval --out [OUT]
`

**1**. `--input-h5ad [INPUT_H5AD]` specifies an h5ad file name (from the AnnData package) containing raw gene expression counts for all cells. Sample-level covariates should be found in the `obs` variable (see [here](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) for info on AnnData objects).   

**2**. `--perturbation-column-name [PERTURBATION_COLUMN_NAME]` specifies which column name in the cell metadata matrix (stored in `.obs` in the h5ad object) corresponds to perturbation information. The values in this column should represent which perturbations (by name) are present in each cell. In the below example, the column `Perturbation` contains this info. 

```
Cell_Barcode  Perturbation  Covariate_1  Covariate_2  ...
Cell_1  Pert_1  0.04  45
Cell_2  Control  0.14  69
Cell_3  Pert_2  0.29  12
Cell_4  Pert_1:Pert_2  0.22  45
...
```

**3**. `--perturbation-delimiter [PERTURBATION_DELIMITER]` (Optional). If a cell contains multiple perturbations, then its value in the perturbation metadata column `.obs[PERTURBATION_COLUMN_NAME]` should correspond to the names of the perturbations separated by some delimiter. In the above example, `Cell_4` contains two perturbations `Pert_1` and `Pert_2` separated by the delimiter `:`. If this flag is set to the delimiter (in this case `:`, but in principle the delimiter can be any string), then multiple perturbations within each cell will be separated according to the delimiter, and effect sizes will be estimated for only individual perturbations. If this flag is not set, then each combination of perturbations is treated as a distinct perturbation, and pertubation effect sizes will be estimated for each distinct perturbation. 

**4**. `--control-perturbation-name [CONTROL_PERTURBATION_NAME]` (Optional). Typically, perturbation effect sizes are computed relative to cells containing a control perturbation, so these control cells should be present in the expression matrix. This flag specifies the perturbation name corresponding to a control guide. The name must be present in `.obs[PERTURBATION_COLUMN_NAME]`. If multiple different names correspond to control guides, then the names can be provided as a list separated by commas, e.g. `ctrl_pert1,ctrl_pert2,ctrl_pert3`. If this flag is not set, then perturbation effect sizes will be computed relative to the mean expression across all cells. 

**5**. `--covariates [COVARIATES]` (Optional). This optional flag specifies a comma-separated list of covariate names to regress out of the expression matrix, e.g. `cov_1,cov_2,cov_3`. These names must be present in the column names of the `obs` variable in the h5ad file. 

**6**. `--compute-pval` (Optional). This optional flag specifies whether or not to compute p-values for all effect sizes by permutation testing. If set, two additional files with suffixes `_pvals.txt` and `_qvals.txt` will be output. 

**7**. `--out [OUT]` specifies the output file prefix, including the directory. The output name for the effect sizes will be the prefix with `_LFCs.txt` appended to it. 
<br /><br />
# Output

FR-Perturb will output a whitespace-delimited table with suffix `_LFCs.txt` containing the log-fold changes of expression relative to expression in cells containing control guides (if `--control-perturbation-name` is set), or relative to mean expression across all cells (if `--control-perturbation-name` is not set). Rows correspond to genes and columns correspond to perturbations. An example is shown below:

```
Pert_1  Pert_2  Pert_3
Gene_1  0.39 -0.9 2.1 
Gene_2  0.32 -0.82 3.1
Gene_3  -0.61 0.24 1.2
```

If `--compute-pval` is set, then FR-Perturb will also output additional two files with suffixes `_pvals.txt` and `_qvals.txt` containing p-values (with the same rows and columns as the effect size file) and q-values (computed using the Benjamini-Hochberg procedure) respectively. 
<br /><br />

# Options for large datasets
If your Perturb-seq dataset is large (i.e. greater than around 250,000 cells) and you're having trouble running FR-Perturb due to memory issues, or if FR-Perturb is taking too long to compute p-values, then the following flags may be useful.

**1**. `--large-dataset`. This flag should be set if memory is an issue when processing the dataset using FR-Perturb with default parameters. With 100G of RAM, default FR-Perturb can handle datasets with ~250K cells. If this flag is set, FR-Perturb can handle datasets with up to millions of cells, though compute time is substantially increased. It uses memory-efficient randomized partial SVD rather than sparse PCA during factorize step of factorize-recover. 

**2**. `--batches [BATCHES]`. If `--large-dataset` is set, then data is processed in batches. This flag specifies how many batches to use. Increasing this number reduces the amount of memory but increases compute time. Default 10. 

**3**. `--temp-out [TEMP_OUT]`. Set this flag if you want to compute p-values in separate batches, which you might want to do if computing p-values in the main script takes too long. These batches can be run in parallel to reduce running time, but you need to implement the parallelization yourself (sorry). If you use a compute cluster then this can be easily done by submitting each batch as a separate job. There are three steps: 
- **Step 1.** Run the main script with the `--temp-out` flag and specify a temporary output prefix (including the full path) to save a temporary file needed to perform permutations. Example:
  - `./run_FR_Perturb.py --input-h5ad [INPUT_H5AD] --input-perturbation-matrix [INPUT_PERTURBATION_MATRIX] --temp-out [TEMP_OUT] --out [OUT]`. 
- **Step 2.** Run multiple instances of the compute_pvalues_batched.py script in parallel. Example:
  - `for BATCH_NUMBER in 1:NUM_BATCHES (parallelize this however you'd like); do ./compute_pvalues_batched.py --input-prefix [TEMP_OUT] --batch-number [BATCH_NUMBER] --perms-per-batch [PERMS_PER_BATCH]`
  - `--input-prefix [TEMP_OUT]`. This argument should be exactly the same as `[TEMP_OUT]` from Step 1.
  - `--batch-number [BATCH_NUMBER]`. This should be some unique number for the given batch of permutations. You can set it to whatever number you like as long as each batch gets a unique number. 
  - `--perms-per-batch [PERMS_PER_BATCH]`. This should be an integer that specifies how many permutations to perform for the given batch. The total number of perturbations will be the total number of times you run the script multipled by the number of permutations per batch. A reasonable option is to run 40 instances of the compute_pvalues_batched.py script with 25 permutations per batch, resulting in 1000 total pertutations. 
- **Step 3.** After all the previous scripts have finished, run the compute_pvalues batched.py script one more time to combine the permutation files and compute p-values:
  - `./compute_pvalues_batched.py --input-prefix [TEMP_OUT] --combine-pval-files --out [OUT]`
  - `--input-prefix [TEMP_OUT]`. This argument should be exactly the same as `[TEMP_OUT]` from Step 1.
  - `--out [OUT]`. Output file prefix for the pvalues, including the directory. This can match `[OUT]` from Step 1, or you can set it to something different if you want. 

To summarize, this is an example sequence of commands you should run to compute p-values in batches:
```
./run_FR_Perturb.py --input-h5ad [INPUT_H5AD] --input-perturbation-matrix [INPUT_PERTURBATION_MATRIX] --temp-out [TEMP_OUT] --out [OUT]

## Submit each instance of the loop as a separate job, or parallelize this however else you'd like
for BATCH_NUMBER in {1..40}; do
  ./compute_pvalues_batched.py --input-prefix [TEMP_OUT] --batch-number [BATCH_NUMBER] --perms-per-batch 25
done

./compute_pvalues_batched.py --input-prefix [TEMP_OUT] --combine-pval-files --out [OUT]
```


# Other options

FR-Perturb has additional options that can be used to tweak data processing, model hyperparameters, p-value calculation, and other things. 

### Gene names

**1**. `--gene-column-name [GENE_COLUMN_NAME]`. By default, FR-Perturb will assume that the index of `.var` of the AnnData object represents gene names when outputting results. If a different column of `.var` corresponds to the desired gene names, then the column name should be specified with this flag. 

### Significance testing parameters 

**1**. `--num-perms [NUM_PERMS]` specifies the number of permutations to perform for permutation testing. Default 1000. 

**2**. `--fit-zero-pval` specifies whether or not to fit a skew-normal distribution the null distribution of effect sizes with p=0 from permutation testing. Setting this flag allows for p-values below 1/num_perms, though it substantially increases compute time. 

**3**. `--multithreaded` specifies whether or not use multithreading when fitting skew-normal distributions. 

### Model hyperparameters

**1**. `--rank [RANK]` specifies the rank of the matrix during the factorize step of factorize-recover. Default 20. We tested this number across a bunch of Perturb-seq datasets and it's sensible.  

**2**. `--spca-alpha [SPCA_ALPHA]` specifies a hyperparameter determining the sparsity of the factor matrix during the factorize step of factorize-recover. Higher value = more sparse. Default 0.1. We tested this number across a bunch of Perturb-seq datasets and it's sensible.  

**3**. `--lasso-alpha [LASSO_ALPHA]` specifies a hyperparameter determining the sparsity of learned effects during the recover step of factorize-recover. Higher value = more sparse. Default 0.0001. We tested this number across a bunch of Perturb-seq datasets and it's sensible.  

### Running other perturbation effect inference methods

**1**. `--elastic-net`. Estimate effect sizes using elastic net as done in Dixit et al. 2016 Cell rather than FR-Perturb. 

**2**. `--elastic-net-genes-per-batch`. Optional argument to split up the inference over batches of downstream genes to reduce memory usage. Input is the number of genes in each batch (larger = fewer batches but more memory).

### Cross-validation

**1**. `--cross-validate`. FR-Perturb provides the option to perform 2-fold cross validation of the perturbation effects. The input dataset is split in half, then perturbation effect sizes are separately estimated in each half and their consistency with each other is assessed. If `--compute-pval` is not set, then correlation and sign consistency will be compute for the largest effects by magnitude. The output will look something like this:
```
        Correlation     Sign_concordance
Top_0.1%_effects        0.74     0.87
Top_1.0%_effects        0.45     0.77
```

If `--compute-pval` is set, then correlation and sign consistency will be computed for the largest effects by magnitude as well as all significant effects. The output will look something like this:
```
        Correlation     Sign_concordance
Top_0.1%_effects        0.74     0.87
Top_1.0%_effects        0.45     0.77
q<0.05_effects  0.83     0.95
q<0.2_effects   0.62      0.82
```

**2**. `--cross-validate-runs [CROSS_VALIDATE_RUNS]`. Number of times to repeat cross-validation. Cross-validation accuracy is averaged across runs. Default 2. 






