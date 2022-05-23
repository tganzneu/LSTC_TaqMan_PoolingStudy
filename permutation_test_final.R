library(dplyr)
# read in data for asymptomatic samples and clean
tf_asymp_R <- read.csv("~/neu/Final Reports/Thermo Pooling/PlosOne/R_analysis/tf_asymp_R.csv")
tf_asymp_R <- tf_asymp_R[-c(1,8:10)]
View(tf_asymp_R)
# get means and medians for each column
a_means <- as.data.frame(t(apply(tf_asymp_R, 2, mean, na.rm=TRUE)))
a_meds <- as.data.frame(t(apply(tf_asymp_R, 2, median, na.rm=TRUE)))
# calculate test statistics: absolute difference between pool and neat 
# for means and medians
stat_means <- c()
stat_meds <- c()
for (i in (1:3)){
  stat_means[i] = ((a_means[i]) - (a_means[i+3]))
  stat_meds[i] = ((a_meds[i]) - (a_meds[i+3]))
}
stat_means <- unlist(stat_means)
stat_meds <- unlist(stat_meds)
a_test_stats <- data.frame(stat_means, stat_meds)   
row.names(a_test_stats) <- c("O","N","S")
a_test_stats
# Pivot table code to set up dataframe for permutation
a_o <- na.omit(data.frame(stack(tf_asymp_R[c(1,4)])))
a_n <- na.omit(data.frame(stack(tf_asymp_R[c(2,5)])))
a_s <- na.omit(data.frame(stack(tf_asymp_R[c(3,6)])))
# Define function to permute
permute <- function(data, permutations){
  set.seed(10)
  permsamp <- matrix(0, nrow=nrow(data), ncol=permutations)
  for (i in 1:permutations){
    permsamp[,i] <- sample(data$values, size=nrow(data), replace=FALSE)
  }
  return(permsamp)
}
# Permute for each gene for number of cycles
a_o_p <- permute(a_o, 1e6)
a_n_p <- permute(a_n, 1e6)
a_s_p <- permute(a_s, 1e6)
# Define function to evaluate test statistics
perm_stat <- function(data, permutations){
  p_stat_mean <- p_stat_med <- rep(0,permutations)
  for (i in 1:permutations){
    p_stat_mean[i] <- (mean(data[1:(nrow(data)/2),i])-
                            mean(data[((nrow(data)/2)+1):(nrow(data)),i]))
    p_stat_med[i] <- (median(data[1:(nrow(data)/2),i])-
                            median(data[((nrow(data)/2)+1):(nrow(data)),i]))
  }
  df <- data.frame(p_stat_mean, p_stat_med)
  return(df)
}
# Generate test statistics from permuted data
a_o_p_t <- perm_stat(a_o_p, 1e6)
a_n_p_t <- perm_stat(a_n_p, 1e6)
a_s_p_t <- perm_stat(a_s_p, 1e6)
# Calculate p-values for both means and medians
a_o_p_t_mn <- mean(a_o_p_t$p_stat_mean >= a_test_stats$stat_means[1])
a_n_p_t_mn <- mean(a_n_p_t$p_stat_mean >= a_test_stats$stat_means[2])
a_s_p_t_mn <- mean(a_s_p_t$p_stat_mean >= a_test_stats$stat_means[3])
a_o_p_t_md <- mean(a_o_p_t$p_stat_med >= a_test_stats$stat_meds[1])
a_n_p_t_md <- mean(a_n_p_t$p_stat_med >= a_test_stats$stat_meds[2])
a_s_p_t_md <- mean(a_s_p_t$p_stat_med >= a_test_stats$stat_meds[3])
a_pvals <- data.frame(c(a_o_p_t_mn,a_n_p_t_mn,a_s_p_t_mn),
                      c(a_o_p_t_md,a_n_p_t_md,a_s_p_t_md))
colnames(a_pvals) <- c('mean','median')
row.names(a_pvals) <- c('O','N','S')
a_pvals

# Using the same code, we can generate p values for symptomatic data
# read in data for symptomatic samples and clean
tf_symp_R <- read.csv("~/neu/Final Reports/Thermo Pooling/PlosOne/R_analysis/tf_symp_R.csv")
tf_symp_R <- tf_symp_R[-c(1,8:10)]
View(tf_symp_R)
# get means and medians for each column
s_means <- as.data.frame(t(apply(tf_symp_R, 2, mean, na.rm=TRUE)))
s_meds <- as.data.frame(t(apply(tf_symp_R, 2, median, na.rm=TRUE)))
# calculate test statistics: absolute difference between pool and neat 
# for means and medians
stat_means <- c()
stat_meds <- c()
for (i in (1:3)){
  stat_means[i] = ((s_means[i]) - (s_means[i+3]))
  stat_meds[i] = ((s_meds[i]) - (s_meds[i+3]))
}
stat_means <- unlist(stat_means)
stat_meds <- unlist(stat_meds)
s_test_stats <- data.frame(stat_means, stat_meds)   
row.names(s_test_stats) <- c("O","N","S")
s_test_stats
# Pivot table code to set up dataframe for permutation
s_o <- na.omit(data.frame(stack(tf_symp_R[c(1,4)])))
s_n <- na.omit(data.frame(stack(tf_symp_R[c(2,5)])))
s_s <- na.omit(data.frame(stack(tf_symp_R[c(3,6)])))
# Permute for each gene for number of cycles
s_o_p <- permute(s_o, 1e6)
s_n_p <- permute(s_n, 1e6)
s_s_p <- permute(s_s, 1e6)
# Generate test statistics from permuted data
s_o_p_t <- perm_stat(s_o_p, 1e6)
s_n_p_t <- perm_stat(s_n_p, 1e6)
s_s_p_t <- perm_stat(s_s_p, 1e6)
# Calculate p-values for both means and medians
s_o_p_t_mn <- mean(s_o_p_t$p_stat_mean >= s_test_stats$stat_means[1])
s_n_p_t_mn <- mean(s_n_p_t$p_stat_mean >= s_test_stats$stat_means[2])
s_s_p_t_mn <- mean(s_s_p_t$p_stat_mean >= s_test_stats$stat_means[3])
s_o_p_t_md <- mean(s_o_p_t$p_stat_med >= s_test_stats$stat_meds[1])
s_n_p_t_md <- mean(s_n_p_t$p_stat_med >= s_test_stats$stat_meds[2])
s_s_p_t_md <- mean(s_s_p_t$p_stat_med >= s_test_stats$stat_meds[3])
s_pvals <- data.frame(c(s_o_p_t_mn,s_n_p_t_mn,s_s_p_t_mn),
                      c(s_o_p_t_md,s_n_p_t_md,s_s_p_t_md))
colnames(s_pvals) <- c('mean','median')
row.names(s_pvals) <- c('O','N','S')
s_pvals
# P-values
a_pvals
s_pvals
# Visualize histograms of permuted sampling distributions and show observed value intersect
par(mfrow=c(2,3))
hist(a_o_p_t$p_stat_mean)
abline(v=a_test_stats$stat_means[1])
hist(a_n_p_t$p_stat_mean)
abline(v=a_test_stats$stat_means[2])
hist(a_s_p_t$p_stat_mean)
abline(v=a_test_stats$stat_means[3])
hist(s_o_p_t$p_stat_mean)
abline(v=a_test_stats$stat_means[1])
hist(s_n_p_t$p_stat_mean)
abline(v=a_test_stats$stat_means[2])
hist(s_s_p_t$p_stat_mean)
abline(v=a_test_stats$stat_means[3])
####
# Confirmatory testing with perm package
library(perm)
permTS(tf_asymp_R$ORF1ab_1.5, tf_asymp_R$ORF1ab_Decon, alternative = "greater")
permTS(as.numeric(na.omit(tf_asymp_R$N.gene_1.5)), as.numeric(na.omit(tf_asymp_R$N.gene_Decon)), alternative = "greater")
permTS(as.numeric(na.omit(tf_asymp_R$S.gene_1.5)), as.numeric(na.omit(tf_asymp_R$S.gene_Decon)), alternative = "greater")
permTS(tf_symp_R$ORF1ab_1.5, tf_symp_R$ORF1ab_Decon, alternative = "greater")
permTS(as.numeric(na.omit(tf_symp_R$N.gene_1.5)), as.numeric(na.omit(tf_symp_R$N.gene_Decon)), alternative = "greater")
permTS(as.numeric(na.omit(tf_symp_R$S.gene_1.5)), as.numeric(na.omit(tf_symp_R$S.gene_Decon)), alternative = "greater")
