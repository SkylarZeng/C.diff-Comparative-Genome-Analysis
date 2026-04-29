######################################################################################################################
# Function: get_pariwise_df_with_location
# This function takes as input:
# 1) snp_dist     - a matrix of inter-sample SNP distances
# 2) site_key     - a 2 column data-frame where the first column is sample ID and the second column is study site 
# 3) ref          - the name of the reference sequence, as a string 
# The output of this function is: 
# A data frame where each row is a unique pair of taxa and 6 columns: 
#   id1       - first taxa id
#   id2       - second taxa id
#   snvs      - count of snvs between first and second taxa
#   site1     - name of location associated with first taxa
#   site2     - name of location associated with second taxa
#   intra     - an indicator of intra or inter location, TRUE == intra
######################################################################################################################

get_pairwise_df_with_location <- function(snp_dist, site_key, ref = NULL) {
  
  # make sure snp_dist is a matrix
  
   if (is.matrix(snp_dist) == FALSE) {
     stop("SNV distance data must be a matrix")
   }
  
  # check to see if reference is still in matrix, and if so, make missing
  
  if (!is.null(ref)) {
    snp_dist[ref,] <- NA
    snp_dist[, ref] <- NA
  }
  
  # delete outliers
  
  snp_dist[snp_dist>1000] <- NA
 
 # delete half of matrix as to not repeat pairs
 
  snp_dist[lower.tri(snp_dist)] <- NA
  
  # create data frame of SNV distances with inter vs intra column
  
  pairs <- melt(snp_dist, value.name = "snvs") %>%  # melt matrix into pairs
    filter(Var1 != Var2, !is.na(snvs)) %>% # filter out intra sample rows and repeat half of matrix
    rename(id1 = Var1, id2 = Var2)
    
  # Add indicator variable about intra or inter site based on site key
  
  pairs$id1 <- as.character(pairs$id1)  # tidy up id vars for joining 
  pairs$id2 <- as.character(pairs$id2)
  
  pairs_with_sites <- left_join(pairs, site_key, by = c("id1" = "id")) %>%
    rename(site1 = site)
  
  pairs_with_sites <- left_join(pairs_with_sites, site_key, by = c("id2" = "id")) %>%
    rename(site2 = site)
  
  pairs_with_sites$intra <- ifelse(pairs_with_sites$site1 == pairs_with_sites$site2, TRUE, FALSE)
  
  return(pairs_with_sites)
  
}

######################################################################################################################
# Function: calc_intra_over_inter_snv_distances
# This function takes as input:
# 1) pairwise_df_with_location     - df in the format of output from the above get_pairwise_df_with_location function
# 2) percentile                    - the percentile you are interested in comparing between the two distributions of SNV distances
# The output of this function is: 
# One summary measure of the percentil SNV distance within a location divided by the precentile SNV distance between locations (a ratio)
######################################################################################################################

calc_intra_over_inter_snv_distances <- function(pairwise_df_with_location, percentile) {
  
  statistic <- quantile(filter(pairwise_df_with_location, group_var == TRUE)$snvs, percentile) / quantile(filter(pairwise_df_with_location, group_var == FALSE)$snvs, percentile)

  return(statistic)

} 

######################################################################################################################
# Function: calc_difference_in_proportion_within_SNV_threshold
# This function takes as input:
# 1) data                         - df with a snv_var (the value of pairwise SNV distances) and a group_var (the variable you want to group by for your comparison)
# 2) snv_var                      - name of the variable that contains the pairwise SNV distance values
# 3) group_var                    - name of the variable that you want to group by
# 4) threshold                    - number of SNVs to use as a threshold
# The output of this function is: 
# The difference in proportion paris less than or equal to the input SNV threshold between the input groups
######################################################################################################################

calc_difference_in_proportion_within_SNV_threshold <- function(data, snv_var, group_var, threshold) {
  x <- mutate(data, within_threshold =  ifelse(data[[snv_var]]<=threshold, TRUE, FALSE))
 # group_by(data[group_var]) %>%
  #summarize(prop_within_threshold = mean(within_threshold))

  abs(mean(x$within_threshold[which(x$ST == 1)]) -  mean(x$within_threshold[which(x$ST == 2)])) * 100

  #summary_df <- mutate(data, within_threshold = ifelse(data[[snv_var]] <= threshold, TRUE, FALSE)) %>%
    #group_by(data[[group_var]]) %>%
    #summarize(prop_within_threshold = mean(within_threshold))
  
  #summary_df[[1,2]] - summary_df[[2,2]]
}

######################################################################################################################
# Function: do_permutations
# This function takes as input:
# 1) pairwise_df_with_location      - df in the format of output from the above get_pairwise_df_with_location function
# 2) var_of_interest_index          - column index of the variable you want to permute
# 3) calculate_statistic             - a function that calculates your statistic of interest
# 4) num                            - number of permutations you want to do 
# 5) percentile                     - the percentile you want to input into your statistic of interest function
# The output of this function is: 
# A dataframe that contains num statistics of interest 
######################################################################################################################

do_permutations <- function(pairwise_df_with_location, var_of_interest_index, num, calculate_statistic, percentile) {
  single_perm <- function(pairwise_df_with_location, var_of_interest_index, calculate_statistic, stat_arg1) {
    perm_df <- pairwise_df_with_location[,1:3]
    perm_df$group_var<- sample(pairwise_df_with_location[,var_of_interest_index])  # do permutation of labels
    perm_statistic <- calculate_statistic(perm_df, percentile)[1]
    return(perm_statistic)
  }

   many_perm_stats <- data.frame(replicate(num, single_perm(pairwise_df_with_location, var_of_interest_index, calculate_statistic, percentile)))
   colnames(many_perm_stats)[1] <- "stat"
   return(many_perm_stats)
}

######################################################################################################################
# Function: do_permutations_snv_threshold
######################################################################################################################

do_permutations_snv_threshold <- function(pairwise_df_with_location, num, calculate_statistic, snv_var, group_var, threshold) {
  single_perm <- function(pairwise_df_with_location, calculate_statistic, snv_var, group_var, threshold) {
    perm_df <- pairwise_df_with_location[,1:3]
    perm_df$ST<- sample(pairwise_df_with_location$ST)  # do permutation of labels
    perm_statistic <- calculate_statistic(perm_df,snv_var, group_var, threshold)[1]
    return(perm_statistic)
  }
  
  many_perm_stats <- data.frame(replicate(num, single_perm(pairwise_df_with_location, calculate_statistic, snv_var, group_var, threshold)))
  colnames(many_perm_stats)[1] <- "stat"
  return(many_perm_stats)
}


######################################################################################################################
# Function: plot_pairwise_snp_dists_hist
# This function takes as input:
# 1) pairwise_df_with_location     - df in the format of output from the above get_pairwise_df_with_location function
# The output of this function is: 
# An overlapping histogram plotting SNV distances by location 
######################################################################################################################

plot_pairwise_snp_dists_hist <- function(pairwise_df_with_location) {
  ggplot(pairwise_df_with_location, aes(x = snvs, stat(density), fill = intra)) +
    geom_histogram(alpha = 0.3, position = "identity", binwidth = 1)
}
