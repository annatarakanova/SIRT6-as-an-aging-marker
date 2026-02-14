#!/usr/bin/env Rscript

########################################################
## Differential expression for human SIRT6 datasets ####
########################################################

library(arrow) # read/write parquet files
library(dplyr)
library(tibble)
library(DESeq2)

#########################
## Paths (arguments) ####
#########################

args <- commandArgs(trailingOnly = TRUE) # reads arguments passed to the script via a command line 

get_arg <- function(x) {
  i <- which(args == x)
  if (length(i) == 0) return(NULL)
  args[i + 1]
}

organism <- get_arg("--organism") # The organism
expr_path <- get_arg("--expr_path") # Count matrix
meta_path <- get_arg("--meta_path") # Meta data
out_dir <- get_arg("--out_dir") # Output dir

# Prevent from mistake
if (is.null(expr_path) || is.null(meta_path) ||
    is.null(out_dir) || is.null(organism)) {
  stop("Usage: --expr_path --meta_path --out_dir --organism")
}

# Create output dir if they are not already existed 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), showWarnings = FALSE)

#########################
###### Load metadata ####
#########################

samples_meta <- read_parquet(
  file.path(meta_path, "samples.parquet")
) %>% column_to_rownames("sample_id")

sample2exp <- read_parquet(
  file.path(meta_path, "samples_to_experiment.parquet")
)

# Drop columns in meta consisting of NaNs only
samples_meta <- samples_meta[, colSums(!is.na(samples_meta)) > 0]

#########################
###### Discover GSEs ####
#########################

# Paths to the experiment files
expr_files <- list.files(expr_path, pattern = "\\.parquet$", full.names = TRUE)
gse_ids <- sub("\\.parquet$", "", basename(expr_files))

#########################
# Filtering function ####
#########################

# Filter lowly expressed genes
# The filtering rule: gene is considered "expressed" in a group if it has at least 10 counts in enough samples of that group.
# Gene is kept if it passes these criteria in both of biological groups.
filter_counts <- function(counts, meta, group_var, min_count = 10) { # group_var - a column in meta that defines the groups
  
  # Start by assuming all genes are kept -> then remove genes that fail the filter
  keep <- rep(TRUE, nrow(counts))
  
  # Loop over biological groups 
  for (group in levels(meta[[group_var]])) { 
    # Samples belonging to this group
    s <- rownames(meta) [meta[[group_var]] == group]
    
    # Minimum number of samples required to express the genes
    n_min <- max(2, floor(length(s) / 2)) # at least 2 or half of the samples in a group (whichever is larger)
    
    # Check expression criterion for this group
    keep_group <- rowSums(
      counts[, s, drop = FALSE] >= min_count
      ) >= n_min
    
    # Conservative rule: gene must pass in both of the groups
    keep <- keep & keep_group
  }
  
  # Return filtered count matrix
  counts[keep, , drop = FALSE]
}

#########################
########## Main loop ####
#########################

message("Running DE analysis for organism: ", organism)

for (i in seq_along(expr_files)) { # One iteration = one GSE
  
  gse <- gse_ids[i]
  message("Processing ", gse)
  
  # Create experiment-specific output directory
  gse_out_dir <- file.path(out_dir, "results", gse)
  dir.create(gse_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    
    # Load count matrix
    counts <- read_parquet(expr_files[i]) %>%
      column_to_rownames("gene_id")
    
    # Select only samples for this GSE
    samples <- sample2exp %>%
      filter(experiment_id == gse) %>%
      pull(sample_id)
    
    # Subsets metadata to match the count matrix
    meta <- samples_meta[samples, , drop = FALSE]
    
    ##############################
    ## Define stratification #####
    ##############################
    
    # Split experiments by confounding factors (treatment/cell type/condition/tissue) -> run separate DE per treatment/cell type/condition/tissue
    
    strata <- list(all = seq_len(nrow(meta))) # default: no splitting
    
    # Candidate splitting variables
    strat_vars <- c("treatment", "cell_type", "condition", "tissue")
    
    # Normalize stratification variables
    for (v in strat_vars) {
      if (v %in% colnames(meta)) {
        meta[[v]] <- trimws(as.character(meta[[v]]))
      }
    }
    
    for (v in strat_vars) {
      if (v %in% colnames(meta)) {
        
        x <- meta[[v]]
        
        # Use this variable for splitting only if:
        # 1) not all NA
        # 2) more than one unique non-NA value
        if (!all(is.na(x)) && length(unique(x[!is.na(x)])) > 1) {
          
          message("Stratifying experiment ", gse, " by ", v)
          
          strata <- split(seq_len(nrow(meta)), x)
          break # use only one splitting variable
        }
      }
    }
    
    
    ############################
    ### Run DE per stratum #####
    ############################
    
    for (s in names(strata)) { # s is the name of stratum (DMSO, RAFi, etc.)
      
      idx <- strata[[s]] # indices of samples belonging to the stratum
      meta_s <- meta[idx, , drop = FALSE] # subset metadata to only these samples
      
      counts_s <- counts[, rownames(meta_s), drop = FALSE] # subset the count matrix to the same samples
      
      # Lists all genotypes present
      genotypes <- unique(meta_s$genotype)
      
      message(gse, " genotypes: ", paste(genotypes, collapse = ", "))
      
      # Require WT as a reference
      if (!"WT" %in% genotypes) next
      
      # Genotypes to compare against WT
      test_genotypes <- setdiff(genotypes, "WT")
      
      ############################
      ### Loop over genotypes ####
      ############################
      
      for (g in test_genotypes) { # Each iteration = one contrast (WT vs KO, WT vs Het, etc.)
        message("Testing contrast: WT vs ", g)
        
        # Sanitize genotype name
        g_factor <- make.names(g)
        
        # Subset for contrast
        meta_g <- meta_s %>%
          dplyr::filter(genotype %in% c("WT", g))
        
        message("Sample counts:")
        print(table(meta_g$genotype))
        
        if (any(table(meta_g$genotype) < 2)) next
        
        # Subset counts to the same samples
        counts_g <- counts_s[, rownames(meta_g), drop = FALSE]
        
        # Define comparison variable
        meta_g$genotype_cmp <- factor(
          make.names(meta_g$genotype),
          levels = c("WT", g_factor) # WT is the reference
        )
        
        # Filter lowly expressed genes
        counts_f <- filter_counts(
          counts = counts_g,
          meta = meta_g,
          group_var = "genotype_cmp"
        )
        
        message("Genes after filtering: ", nrow(counts_f))
        
        # Skip contrasts with too few genes after filtering
        if (nrow(counts_f) < 100) next
        
        #######################################
        # Build design formula dynamically ####
        #######################################
        
        # Candidate covariates to include if available
        candidate_covariates <- c("sex", "age", "strain")
        
        usable_covariates <- c()
        
        for (v in candidate_covariates) {
          if (v %in% colnames(meta_g)) {
            x <- meta_g[[v]]
            
            # Use covariates only if:
            # 1) not all NA
            # 2) more than one unique non-NA value
            if (!all(is.na(x)) && length(unique(x[!is.na(x)])) > 1) {
              usable_covariates <- c(usable_covariates, v)
            }
          }
        }
        
        # Sanitize covariate levels (make them R-safe)
        for (v in usable_covariates) {
          meta_g[[v]] <- factor(make.names(meta_g[[v]]))
        }
        
        # Final design terms:
        # covariates first, genotype last (as it's effect of interest)
        design_terms <- c(usable_covariates, "genotype_cmp")
        
        design_formula <- as.formula(
          paste("~", paste(design_terms, collapse = " + "))
        )
        
        message("Design formula: ", deparse(design_formula))
        
        # Build DESeq object
        
        dds <- DESeqDataSetFromMatrix(
          countData = round(counts_f),
          colData = meta_g,
          design = design_formula
        )
        
        # Run DESeq (fits the negative binomial model)
        dds <- DESeq(dds, quiet = TRUE)
        
        # Extract DE results for g vs WT
        res <- results(
          dds,
          contrast = c("genotype_cmp", g_factor, "WT")
        )
        
        ############################
        ##### Collect results ######
        ############################
        
        # Convert DESeq2 results into a tidy data frame
        out <- as.data.frame(res) %>%
          rownames_to_column("gene_id")
        
        # Add a column indicating significant genes
        out <- out %>%
          mutate(
            significant = !is.na(padj) & padj < 0.05
          )
        
        # Annotate results
        out$experiment_id <- gse
        out$stratum <- s
        out$organism <- organism
        out$contrast <- paste0(g_factor, "_vs_WT")
        
        ############################
        ##### Write the output #####
        ############################
        
        # Save files in a safer format (replace anything that is not letters, numbers, ., _, - with _)
        s_safe <- gsub("[^A-Za-z0-9._-]", "_", s)
        g_file <- gsub("[^A-Za-z0-9._-]", "_", g)
        
        # Save DE results in parquet format
        write_parquet(
          out,
          file.path(
            gse_out_dir,
            paste0(gse, "_", s_safe, "_", g_file, "_vs_WT_deseq2.parquet")
          )
        )
      }
    }
    
  # Error handling (if any experiment failed)
  }, error = function(e) {
    cat(
      paste(Sys.time(), gse, e$message, "\n"),
      file = file.path(out_dir, "logs", "errors.log"),
      append = TRUE
    )
  })
}

message(organism, " DE analysis complete.")
