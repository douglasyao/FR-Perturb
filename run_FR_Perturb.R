#!/usr/bin/env Rscript

Sys.setenv(KMP_WARNINGS = 'off')
suppressMessages(library(optparse))
suppressMessages(library(Seurat))
suppressMessages(library(spams))
suppressMessages(library(Matrix))
suppressMessages(library(parallel))
suppressMessages(library(pbapply))

option_list <- list( 
### Required
    make_option("--input_seurat", type ="character",
        help="Input file name containing a Seurat object (stored as an R object) containing the following information:
(1) Raw gene x cell expression count matrix under an assay named 'RNA'
(2) Binary cell x perturbation matrix (with 1 indicating that a cell containing a perturbation, otherwise 0) under an assay named 'perturbations'
(3) Meta-data for each cell"),
    make_option("--control_perturbation_name", type="character",
        help="Comma-separated list of perturbation names that represent control perturbations (names must be present in rows of 'perturbation' assay in Seurat object)" ),
    make_option("--out", type="character",
        help="Output prefix (including directory) for effect sizes (and optionally p-values)"),

### Optional 
    make_option("--compute_pval", action="store_true", 
        help="Whether or not to compute p-values for all effect size estimates by permutation testing"),            
    make_option("--rank", type="integer", default=20, 
        help="Hyperparameter determining the rank of the matrix during the factorize step"),
    make_option("--lambda1", type="numeric", default=0.1, 
        help="Hyperparameter determining the sparsity of the factor matrix during the factorize step of the method. Higher value = more sparse."),
    make_option("--lambda2", type="numeric", default=10, 
        help="Hyperparameter determining the sparsity of learned effects during the recover step of the method. Higher value = more sparse."),
    make_option("--covariates", type="character", 
        help="Comma-separated list of covariate names to regress out of the expression matrix (names must match the column names in the meta-data of the Seurat object)"),
    make_option("--guide_overloaded", action="store_true",
        help="Runs the version of FR-Perturb that assumes data is generated from guide-overloading"),
    make_option("--droplet_overloaded", action="store_true", 
        help="Runs the version of FR-Perturb that assumes data is generated from droplet-overloading"),
    make_option("--num_perms", type="integer", default=10000,
        help="Number of permutations when doing permutation testing"),
    make_option("--fit_zero_pval", action="store_true",
        help="Compute p-values by fitting skew-t distribution to null distribution (allows for p-values below 1/num_perms, but significantly increases compute time)"),
    make_option("--output_factor_matrices", action="store_true",
        help="Whether or not to output the latent gene expression factor matrices in addition to the full effect sizes matrix")
)

regress_covariates <- function(exp.mat, cov.mat, intercept = TRUE, add.back.intercept = FALSE) {
#     Regress covariates from expression matrix 
    
#     Input 
#     --------
#     exp.mat: cell x gene expression matrix
#     cov.mat: cell x covariate matrix
#     intercept: Whether to include intercept when regressing covariates
#     add.back.intercept: Whether to add back intercept to values after regressing out covariates (used if we don't want/expect expression matrix to be centered)
    
#     Returns
#     --------
#     Expression matrix with covariates regressed out
    
    cov.mat <- t(t(cov.mat) - colMeans(cov.mat)) # center covariates
    if (intercept) cov.mat <- cbind(Intercept = rep(1, nrow(cov.mat)), cov.mat)
    cov.mat <- as.matrix(cov.mat)
    
    fit <- .lm.fit(cov.mat, exp.mat)
    
    if (add.back.intercept) exp.mat <- t(t(fit$residuals) + fit$coefficients[1,]) # add back intercept values to points 
    else exp.mat <- fit$residuals
    
    return(exp.mat)
}

scale_effs <- function(B, mean.exp, downsample.num = 25000, log.exp.baseline = 2, plot = FALSE) {
#     Scale effect sizes based on expression. Fits relationship between log expression and log(abs(LFC)).
    
#     Input
#     --------------
#     B: Gene x perturbation effect size matrix
#     mean.exp: Named vector indicating log mean expression of all genes
#     downsample.num: How many effects to downsample to when fitting curve
#     log.exp.baseline: Log expression to scale mean to
#     plot: Whether to plot effects and curve
    
#     Return
#     -------------
#     List of two items:
#         1. Scaled B matrix
#         2. Per-gene scaling factors
    
    common.genes <- intersect(rownames(B), names(mean.exp))
    B <- B[common.genes,]
    mean.exp <- mean.exp[common.genes]
    rand.idx <- sample.int(length(B), downsample.num)
    rand.idx.array <- arrayInd(rand.idx, .dim = dim(B))
    to.plot <- data.frame(x = mean.exp[rand.idx.array[,1]], y = log(abs(B[rand.idx])))
    to.plot$y[is.infinite(to.plot$y)] <- NA
    to.plot <- na.omit(to.plot)
    fit <- loess(y ~ x, data = to.plot)
    to.plot$line <- fit$fitted
    baseline <- to.plot[which.min(abs(to.plot$x - log.exp.baseline)), 'line']
    scale.factors <- exp(predict(fit, mean.exp) - baseline)
    B <- B / scale.factors
    if (plot) {
        to.plot <- to.plot[order(to.plot$x),]
        plot(to.plot$x, to.plot$y); lines(to.plot$x, to.plot$line, col = 'red')
    }
    return(list(B = B, scale.factors = scale.factors)) 
}

fit_skew_t <- function(t_nulls, t_star, side = 'both', max.time = 20) {
#     Fit skew t distribution to empirical null distribution and computes pvalues. Code taken from SCEPTRE (Katsevich et al. 2021 Genome Biol)
    
#     Input
#     -------------------
#     t_nulls: Vectors corresponding to values of null distribution
#     t_star: Observed statistic
#     side: Which sided pvalue to compute ('left', 'right', or 'both')
    
#     Returns
#     ---------------------
#     List of three items
#         1. skew_t_fit_success: Whether the fitting worked
#         2. out_p: Fitted pvalue. Returns NA if fitting didn't work .
#         3. skew_t_mle: Parameters of fitted skew t distribution. 
    
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(withTimeout(sn::selm(t_nulls ~ 1, family = "ST"), timeout = max.time), error = function(e) return(NA)) # run for a maximum of 20 sec
  # skew_t_fit <- R.utils::withTimeout(sn::selm(t_nulls ~ 1, family = "ST"), timeout = max.time) # run for a maximum of 20 sec
  if (class(skew_t_fit) == "selm") { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- switch(side,
                               'left' = pmax(.Machine$double.eps, sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)), # then compute the skew t-based p-value. pst(x = t_star, dp = dp)
                               'right' = pmax(.Machine$double.eps, 1 - sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
                               'both' = pmax(.Machine$double.eps, sn::pst(x = -abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
                                               (1 - sn::pst(x = abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
      )
    }
  }
  # check if the skew-t fit worked
  skew_t_fit_success <- !is.na(p_value_skew_t)
  if (skew_t_fit_success) {
    out_p <- p_value_skew_t
    skew_t_mle <- dp
  } else {
#     out_p <- switch(side,
#                     'left' = mean(c(-Inf, t_nulls) <= t_star),
#                     'right' = mean(c(Inf, t_nulls) >= t_star),
#                     'both' = mean(c(Inf, abs(t_nulls)) >= abs(t_star)))
    out_p <- NA
    skew_t_mle <- c(xi = NA, omega = NA, alpha = NA, nu = NA)
  }
  return(list(skew_t_fit_success = skew_t_fit_success, out_p = out_p, skew_t_mle = skew_t_mle))
}

args <- parse_args(OptionParser(option_list=option_list))
# args <- parse_args(OptionParser(option_list=option_list), args = c('--input_seurat=Simulated_seurat.rds',
#                                                                    '--control_perturbation_name=non-targeting',
#                                                                    '--out=temp/out',
#                                                                    '--compute_pval',
#                                                                    '--temp_dir=/n/scratch3/users/d/dwy6/temp/'))

log.file <- paste0(args$out, '.log')
system(paste0('rm ', log.file))
sink(log.file, split = TRUE)
cat('*************************************************\nFactorize-Recover for Perturb-seq analysis (FR-Perturb)\nDouglas Yao 2022\n*************************************************\n\n')
cat(paste('Call:\n./run_FR_Perturb.R', paste(commandArgs(trailingOnly=TRUE), collapse = ' '), '\n\n'))

if (!is.null(args$guide_overloaded) & !is.null(args$droplet_overloaded)) {
    stop('Only one of --guide_overloaded and --droplet-overloading should be set')
} else if (!is.null(args$droplet_overloaded)) {
    overload.type <- 'droplet'
} else {
    overload.type <- 'guide'
}

# load data
start.time <- Sys.time()
cat(paste0('Analysis started at ', start.time, '\n'))
cat('Loading Seurat data...  ')
dat <- readRDS(args$input_seurat)
cat('Done\n')

# center rows of p.mat
p.mat <- as.matrix(t(dat@assays$perturbations@data))
if (overload.type == 'droplet') p.mat <- p.mat / rowSums(p.mat)

# extract covariates from Seurat meta data
if (!is.null(args$covariates)) {
    cov.names <- strsplit(args$covariates, split = ',')[[1]]
    cov.mat <- dat@meta.data[,cov.names,drop=F]
    unique.count <- sapply(cov.mat, FUN = function(x) length(unique(x)))
    cov.mat <- cov.mat[,!unique.count == 1] # remove covariates that are the same across all samples
    cov.mat[,sapply(cov.mat, is.character)] <- sapply(cov.mat[,sapply(cov.mat, is.character)], factor) # convert characters to factors
    cov.mat <- model.matrix(~ ., data=cov.mat, contrasts.arg=lapply(cov.mat[,sapply(cov.mat, is.factor),drop=F], contrasts, contrasts = TRUE)) # convert factors to model matrix
    cov.mat <- cov.mat[,-1,drop=F]
}

# regress out covariates and center expression matrix based on control expression
cat('Regressing out covariates and centering expression matrix...  ')                           
dat <- NormalizeData(dat, normalization.method = 'RC', verbose = F)
logmeanexp <- log(rowMeans(dat@assays$RNA@data))
dat <- NormalizeData(dat, verbose = F)
if (!is.null(args$covariates)) {
    exp.mat <- as.matrix(t(dat@assays$RNA@data))
    exp.mat <- regress_covariates(exp.mat, cov.mat)
} else {
    exp.mat <- as.matrix(dat@assays$RNA@data)
    exp.mat <- t(exp.mat - rowMeans(exp.mat))
}
n.guides <- colSums(dat[['perturbations']])
nt.names <- strsplit(args$control_perturbation_name, split = ',')[[1]]                           
ctrl.idx <- which((colSums(dat@assays$perturbations@data[rownames(dat[['perturbations']]) %in% nt.names,,drop=F]) == n.guides) & n.guides > 0) # cells containing only control guides
ctrl.exp <- colMeans(exp.mat[ctrl.idx,])
exp.mat <- exp.mat - rep(ctrl.exp, each = nrow(exp.mat)) # center expression matrix based on control expression
cat('Done\n')
    
# Factorize
cat('Factorizing expression matrix...  ')                           
exp.mat <- exp.mat[,!colSums(exp.mat) == 0]
keep.cells <- which(rowSums(p.mat) > 0 & !is.na(rowSums(p.mat))) # factorize with entire expression matrix, recover with only cells w guides
p.mat <- p.mat[keep.cells,]      
exp.mat <- t(exp.mat)
W <- as.matrix(spams.trainDL(exp.mat, K = args$rank, lambda1 = args$lambda1, iter = 50, verbose = FALSE))
U.tilde <- as.matrix(spams.lasso(exp.mat, D = W, lambda1 = args$lambda1, verbose = FALSE))
U.tilde <- U.tilde[,keep.cells]
U.tilde <- t(U.tilde)
W <- t(W)
cat('Done\n')
                           
# Recover
cat('Regressing left factor matrix on perturbation design matrix...  ') 
U <- as.matrix(spams.lasso(U.tilde, D = p.mat, lambda1 = args$lambda2, verbose = FALSE))
B <- t(U %*% W)
colnames(B) <- colnames(p.mat)
rownames(B) <- rownames(exp.mat)    
cat('Done\n')   

# Compute pvalues
if (!is.null(args$compute_pval)) {
    temp.B <- t(B)
    
    if (is.null(args$fit_zero_pval)) {
        cat(sprintf('Computing p-values by permutation testing (%d total permutations)... \n', args$num_perms))
        current.pct <- 0.05
        pvals <- matrix(0L, nrow = nrow(temp.B), ncol = ncol(temp.B))
        pb <- timerProgressBar()
        for (i in 1:args$num_perms) {
            if (i > current.pct * args$num_perms + 1) {
                setTimerProgressBar(pb, current.pct)
                current.pct <- current.pct + 0.05
            }
            p.mat.perm <- p.mat[sample.int(nrow(p.mat)),]
            U.perm <- as.matrix(spams.lasso(U.tilde, D = p.mat.perm, lambda1 = args$lambda2, verbose = FALSE))
            B.perm <- U.perm %*% W
            temp.indices <- temp.B < B.perm
            pvals[temp.indices] <- pvals[temp.indices] + 1 
        }
        setTimerProgressBar(pb, 1)
        pvals <- (pvals + 1) / (args$num_perms + 1)
        pvals[pvals > 0.5] <- 1 - pvals[pvals > 0.5]
        pvals <- t(2 * pvals)
        rownames(pvals) <- rownames(B)
        colnames(pvals) <- colnames(B)
    } else {
        # fit skew t to null distribution
        args$num_perms <- 500
        cat(sprintf('Computing p-values by permutation testing (%d total permutations)... \n', args$num_perms))
        current.pct <- 0.05

        B.perms <- matrix(0L, nrow = prod(dim(B)), ncol = args$num_perms)
        pb <- timerProgressBar(style = 3, width = 50, char = '+')

        for (i in 1:args$num_perms) {
             if (i > current.pct * args$num_perms + 1) {
                setTimerProgressBar(pb, current.pct)
                current.pct <- current.pct + 0.05
            }
            p.mat.perm <- p.mat[sample.int(nrow(p.mat)),]
            U.perm <- as.matrix(spams.lasso(U.tilde, D = p.mat.perm, lambda1 = args$lambda2, verbose = FALSE))
            B.perms[,i] <- c(U.perm %*% W)
        }
        setTimerProgressBar(pb, 1)
        cat('\n')
        pvals <- rowSums(B.perms < c(temp.B)) / ncol(B.perms)                     
        pvals[pvals > 0.5] <- 1 - pvals[pvals > 0.5]
        pvals <- 2 * pvals

        cat(sprintf('Fitting skew-t distribution to effects with p = 0 (%d total effects)... \n', sum(pvals == 0)))
        zero.indices <- which(pvals == 0)
        B.perms.zero <- B.perms[zero.indices,]
        B.perms.zero <- cbind(temp.B[zero.indices], B.perms.zero)

        n.cores <- detectCores()
        cl <- makeCluster(n.cores)
        clusterExport(cl = cl, 'fit_skew_t')
        invisible(clusterEvalQ(cl = cl, expr = {
            library(R.utils)
            library(sn)
        }))

        pboptions(type = "timer", style = 3)
        new.ps <- pbapply(B.perms.zero, MARGIN = 1, FUN = function(x) fit_skew_t(x[2:length(x)], x[1], side = 'both')$out_p, cl = cl)  
        stopCluster(cl)               
        pvals[zero.indices] <- new.ps
        pvals <- t(matrix(pvals, nrow = ncol(B), ncol = nrow(B)))
        rownames(pvals) <- rownames(B)
        colnames(pvals) <- colnames(B)                  
    }
    pvals <- signif(pvals, digits = 3)
}

# Output files
cat('Outputting results...  ')
B.scaled <- scale_effs(B, logmeanexp)$B
B.scaled <- round(B.scaled, digits = 3)
write.table(B.scaled, file = paste0(args$out, '_effects.txt'), quote = F, sep = '\t')
if (!is.null(args$compute_pval)) write.table(pvals, file = paste0(args$out, '_pvalues.txt'), quote = F, sep = '\t')                         
cat('Done\n')                           

end.time <- Sys.time()                           
cat(paste0('Analysis complete at ', end.time, '\n'))
time.elapsed <- end.time - start.time                        
cat(paste('Total time elapsed:', round(time.elapsed, digits = 2), units(time.elapsed)), '\n')
 