#! Rscript
# epitope_calls_analysis.R

# This script will take pepmeld processed lseq matrix
# 1.) expand to probe-level matrix
# 2.) apply smoothing
# 3.) average replicates
# 4.) create sd-based, pseudo-pval matrix and call epitopes

###############
## Functions ##
###############
source("functions.R")

###########
## Setup ##
###########
# Load libraries
library(dplyr);
library(Matrix);

# load probe_meta for probe_seq <-> probe_id conversion
probe_meta <- readRDS("probe_meta.RDS")

#########
## IgG ##
#########
IgG_sample_meta <- readRDS("IgG_samples.RDS")
IgG_lseq <- read.delim("IgG.lseq.tsv", row.names = 1, check.names = FALSE)
IgG_lseq$PROBE_SEQUENCE <- rownames(IgG_lseq)
IgG_lprobe <- dplyr::left_join(IgG_lseq, dplyr::select(probe_meta, PROBE_SEQUENCE, PROBE_ID), 
                             by = "PROBE_SEQUENCE")

IgG_lprobe_qn <- as.matrix(dplyr::select(IgG_lprobe, -PROBE_SEQUENCE, -PROBE_ID))
rownames(IgG_lprobe_qn) <- IgG_lprobe$PROBE_ID

probe_df <- dplyr::filter(probe_meta, PROBE_ID %in% rownames(IgG_lprobe_qn)) %>%
              dplyr::select(PROBE_ID, SEQ_ID, POSITION)
IgG_lprobe_smooth <- smoothProbeMat(IgG_lprobe_qn,
                                    probe_meta = probe_df, 
                                    dist_method = "aa", 
                                    w = 2, weighted=FALSE)

# average and collapse technical replicates
non_rep_samples <- filter(IgG_sample_meta, is.na(Replicate)) %>% pull(SampleName)
rep_samples <- filter(IgG_sample_meta, !is.na(Replicate)) %>% group_by(Pair) %>% tidyr::nest()
rep_mean <- purrr::map(rep_samples$data, function(df, lprobe_mat) {
                return(rowMeans(lprobe_mat[,df$SampleName], na.rm = TRUE))
            }, lprobe_mat = IgG_lprobe_smooth)
rep_samples <- tidyr::unnest(rep_samples, cols = data) %>% arrange(Pair, Replicate) %>% filter(Replicate == "A")
names(rep_mean) <- rep_samples$SampleName
rep_mean <- as.matrix(bind_cols(rep_mean))
IgG_lprobe_final <- cbind(IgG_lprobe_smooth[,non_rep_samples], rep_mean)

# Update IgG_sample_meta after averaging replicated sample data
IgG_sample_meta <- filter(IgG_sample_meta, SampleName %in% colnames(IgG_lprobe_final))

# calculate probe z-score
IgG_global_mean <- 0.736643108827108
IgG_global_sd <- 	1.03776496459744
IgG_lprobe_zscore <- (IgG_lprobe_final - IgG_global_mean)/IgG_global_sd

# convert z-score to pseudopval for epitope calls
IgG_lprobe_pval <- IgG_lprobe_zscore
IgG_lprobe_pval[IgG_lprobe_pval[,] < 3] <- NA
IgG_lprobe_pval[!is.na(IgG_lprobe_pval)] <- 0.0001
IgG_lprobe_pval[is.na(IgG_lprobe_pval)] <- 1

IgG_epitopes <- makeEpitopeCalls(IgG_lprobe_pval)

# drop control samples before reporting
IgG_sample_meta <- filter(IgG_sample_meta, is.na(Control))

# write epitope calls to xlsx file
results <- purrr::map(IgG_sample_meta$SampleName, function(sample, epitopes) {
                        df <- epitopes[[sample]]
                        df <- select(df, Sample, Epitope_ID, Protein, Start, Stop)
                        colnames(df)[4:5] <- c("Start_probe", "Stop_probe")
                        rownames(df) <- NULL
                        return(df)
                      }, epitopes = IgG_epitopes$sample_epitopes)
results <- bind_rows(results)

fid <- "plasma_IgG_sample_epitopes.xlsx"
openxlsx::write.xlsx(results, file = fid, overwrite = TRUE, sheetName = "epitopes")

#########
## IgA ##
#########
IgA_sample_meta <- readRDS("IgA_samples.RDS")
IgA_lseq <- read.delim("IgA.lseq.tsv", row.names = 1, check.names = FALSE)
IgA_lseq$PROBE_SEQUENCE <- rownames(IgA_lseq)
IgA_lprobe <- dplyr::left_join(IgA_lseq, dplyr::select(probe_meta, PROBE_SEQUENCE, PROBE_ID), 
                               by = "PROBE_SEQUENCE")

IgA_lprobe_qn <- as.matrix(dplyr::select(IgA_lprobe, -PROBE_SEQUENCE, -PROBE_ID))
rownames(IgA_lprobe_qn) <- IgA_lprobe$PROBE_ID

probe_df <- dplyr::filter(probe_meta, PROBE_ID %in% rownames(IgA_lprobe_qn)) %>%
  dplyr::select(PROBE_ID, SEQ_ID, POSITION)
IgA_lprobe_smooth <- smoothProbeMat(IgA_lprobe_qn,
                                    probe_meta = probe_df, 
                                    dist_method = "aa", 
                                    w = 2, weighted=FALSE)

# average and collapse technical replicates
non_rep_samples <- filter(IgA_sample_meta, is.na(Replicate)) %>% pull(SampleName)
rep_samples <- filter(IgA_sample_meta, !is.na(Replicate)) %>% group_by(Pair) %>% tidyr::nest()
rep_mean <- purrr::map(rep_samples$data, function(df, lprobe_mat) {
  return(rowMeans(lprobe_mat[,df$SampleName], na.rm = TRUE))
}, lprobe_mat = IgA_lprobe_smooth)
rep_samples <- tidyr::unnest(rep_samples, cols = data) %>% arrange(Pair, Replicate) %>% filter(Replicate == "A")
names(rep_mean) <- rep_samples$SampleName
rep_mean <- as.matrix(bind_cols(rep_mean))
IgA_lprobe_final <- cbind(IgA_lprobe_smooth[,non_rep_samples], rep_mean)

# Update IgA_sample_meta after averaging replicated sample data
IgA_sample_meta <- filter(IgA_sample_meta, SampleName %in% colnames(IgA_lprobe_final))

# calculate probe z-score
IgA_global_mean <- 0.688089421112556
IgA_global_sd <- 		0.92455735340697
IgA_lprobe_zscore <- (IgA_lprobe_final - IgA_global_mean)/IgA_global_sd

# convert z-score to pseudopval for epitope calls
IgA_lprobe_pval <- IgA_lprobe_zscore
IgA_lprobe_pval[IgA_lprobe_pval[,] < 2.5] <- NA
IgA_lprobe_pval[!is.na(IgA_lprobe_pval)] <- 0.0001
IgA_lprobe_pval[is.na(IgA_lprobe_pval)] <- 1

IgA_epitopes <- makeEpitopeCalls(IgA_lprobe_pval)

# drop control samples before reporting
IgA_sample_meta <- filter(IgA_sample_meta, is.na(Control))

# write epitope calls to xlsx file
results <- purrr::map(IgA_sample_meta$SampleName, function(sample, epitopes) {
  df <- epitopes[[sample]]
  df <- select(df, Sample, Epitope_ID, Protein, Start, Stop)
  colnames(df)[4:5] <- c("Start_probe", "Stop_probe")
  rownames(df) <- NULL
  return(df)
}, epitopes = IgA_epitopes$sample_epitopes)
results <- bind_rows(results)

fid <- "breast_milk_IgA_sample_epitopes.xlsx"
openxlsx::write.xlsx(results, file = fid, overwrite = TRUE, sheetName = "epitopes")

# END

