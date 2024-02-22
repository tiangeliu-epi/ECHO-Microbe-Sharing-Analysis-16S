######      To Create Functions for the ECHO OIF EC0556A Data Analysis      #####
######      Author: Tiange Liu
######      Last updated on 12/1/2024


####=============  Load libraries   ============####
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)
library(readr)
library(tidyverse)
library(stringr)
library(table1)
library(nlme)
library(FEAST)
library(ggVennDiagram)
library(rstatix) # to add significance in ggplot
library(digest)
library(vegan) # to calculate beta diversity matrix


####=============  Below is a list of colors used in figures  =============####
color1 <- c("#BBBBBB", "#E69F00", "#0072B2") #grey, yellow, blue
color1_2 <- c("#BBBBBB", "#0072B2")

color2 <- c("#90C987", "#E69F00", "#0072B2") #green, yellow, blue

color3 <- c("#BBBBBB", "#90C987", "#E69F00", "#0072B2")
color3_2 <- c("white", "#BBBBBB", "#90C987", "#E69F00", "#0072B2")


color4 <- c("#0072B2", "#CC79A7")
color5 <- c("#009E73", "#E69F00")

color6 <- c("#999999", "#CC79A7")

color7 <- c("white","#006A8D", "#B27F10") # white, dark blue, dark yellow


####==============   Write a function to plot split violin plot   =====================####
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}




####=============   Summarize dissimilarity matrix for pairs of mother-own child, and mother-unrelated child  ====================####
# For each dissimilarity matrix, create 4 datasets (stool vs stool) or 3 datasets (stool vs vagina) for each combination of pairs; 
# and then summarize distance matrix       

# Pairs are infant/own mother, infant/unrelated mother, mother/mother, infant/infant (only for stool vs stool)   

cre_dta <- function(matrix, id1, id2){
  dta <- filter(matrix, row.names(matrix) %in% as.list(id1$sample.id))
  names <- row.names(dta)
  
  dta <- as.data.frame(dta[, colnames(dta) %in% as.list(id2$sample.id)])

  row.names(dta) <- names
  dta <- dta[order(row.names(dta)), order(colnames(dta))]
  
  return(dta)
}

sum_diss <- function(betamatrix, infantid, motherid){
  ## infant & mother dataset
  s_i_m <- as.data.frame(cre_dta(betamatrix, infantid, motherid))
  
  
  # Find matching mothers for each infant based on characters 7 to 15
  matching_mother_indices <- sapply(rownames(s_i_m), function(infant_name) {
    infant_substr <- substr(infant_name, 7, 15)  # Extract characters 7 to 15 from the infant name
    matching_mother_index <- match(infant_substr, substr(colnames(s_i_m), 7, 15))
    return(matching_mother_index)
  })
  
  
  ## sum infant vs own mother
  sum1 <- as.data.frame(s_i_m[cbind(1:nrow(s_i_m), matching_mother_indices)])
  colnames(sum1) <- c("values")
  sum1$comparison <- "Child and biological mother"
  
  ## sum infant vs unrelated (random) mother
  for (i in 1:length(matching_mother_indices)) {
    row_index <- i
    col_index <- matching_mother_indices[i]
    s_i_m[row_index, col_index] <- NA
  }
  
  # randomly sample a column index for each row
  num_rows <- nrow(s_i_m)
  
  # initialize a vector to store the randomly selected values
  sum2 <- numeric(num_rows)
  
  # loop through rows and randomly select values, replacing NA with new random values
  set.seed(42)
  for (i in 1:num_rows) {
    repeat {
      random_column_index <- sample(1:ncol(s_i_m), 1)
      row_value <- s_i_m[i, random_column_index]
      if (!is.na(row_value)) {
        sum2[i] <- row_value
        break  # Exit the loop if a non-NA value is found
      }
    }
  }  
  sum2 <- as.data.frame(sum2)
  colnames(sum2) <- c("values")
  sum2$comparison <- "Child and random mother"
  
  
  ## make a summary dataset
  sum_dta <- rbind(sum1, sum2)
  name <- deparse(substitute(betamatrix))
  sum_dta$beta <- ifelse(name=="bc" | name=="bc_2", "Bray-Curtis",
                         ifelse(name=="jc" | name=="jc_2", "Jaccard",
                                ifelse(name=="uu", "Unweighted UniFrac", "Weighted UniFrac")))
  
  return(sum_dta)
}

combine_sum_diss <- function(bc, jc, uu, wu, i_id, m_id){
  sum1 <- sum_diss(bc, i_id, m_id)
  sum2 <- sum_diss(jc, i_id, m_id)
  sum3 <- sum_diss(uu, i_id, m_id)
  sum4 <- sum_diss(wu, i_id, m_id)
  
  sum <- rbind(sum1, sum2, sum3, sum4)
  sum$comparison <- factor(sum$comparison, levels=c("Child and biological mother", "Child and random mother"))
  sum$beta <- factor(sum$beta, levels=c("Jaccard", "Bray-Curtis", "Unweighted UniFrac", "Weighted UniFrac"))
  
  return(sum)
}



####=============  Extract ASV taxonomy and modify name   ==================###
# Input asv table needs to have rows as sample with sample.id as rownames,
  # and columns as taxa with taxa id as colnames. 
# Unnecessary rows need to be removed.

fun_asv_name <- function(asv_table){
  asv_name <- asv_table[c(nrow(asv_table)), ]
  asv_name[1, 1:5]
  
  # convert rows as ASV and columns as names
  asv_name <- as.data.frame(t(asv_name))
  
  # create a separate column to store ASV id
  asv_name$asv <- row.names(asv_name)
  
  # extract names by taxonomy
  junk <- strsplit((asv_name$Taxon), "[;]")
  junk <- as.data.frame(t(sapply(junk, "[", i = 1:max(sapply(junk, length)))))
  
  asv_name$domain <- sub("d__", "", junk$V1)
  asv_name$phylum <- sub("p__", "", junk$V2)
  asv_name$class <- sub("c__", "", junk$V3)
  asv_name$order <- sub("o__", "", junk$V4)
  asv_name$family <- sub("f__", "", junk$V5)
  asv_name$genus <- sub("g__", "", junk$V6)
  asv_name$genus <- sub(" ", "", asv_name$genus)
  asv_name$species <- sub("s__", "", junk$V7)
  asv_name$species <- sub(" ", "", asv_name$species)
  
  asv_name$species <- ifelse(asv_name$genus=="[Bacteroides]_pectinophilus_group", "pectinophilus", asv_name$species)
  asv_name$genus <- ifelse(asv_name$genus=="[Bacteroides]_pectinophilus_group", "Bacteroides", asv_name$genus)
  
  # modify ASV names 
  asv_name$genus <- ifelse(is.na(asv_name$genus) | str_detect(asv_name$genus, "uncultured")==TRUE, "genus unknown", asv_name$genus)
  asv_name$species <- ifelse(is.na(asv_name$species) | str_detect(asv_name$species, "uncultured")==TRUE, "sp. unknown", asv_name$species)
  
  asv_name <- asv_name %>%
    group_by(phylum, class, order, family, genus, species) %>%
    mutate(suffix=1:n())
  
  asv_name$name <- ifelse(asv_name$species!="sp. unknown", paste0(asv_name$genus, " ", asv_name$species, " ", asv_name$suffix),
                          paste0(asv_name$genus, " ", asv_name$suffix))
    # name: genus suffix, or genus species suffix 
  
  asv_name$name <- str_replace_all(asv_name$name, "_", " ")
  
  # generate md5 hash name
  asv_name$md5 <- sapply(asv_name$asv, digest, algo="md5")
  asv_name$name_md5 <- ifelse(asv_name$species!="sp. unknown", paste0(asv_name$genus, " ", asv_name$species, "-", substr(asv_name$md5, 1, 4)),
                              paste0(asv_name$genus, "-", substr(asv_name$md5, 1, 4)))
  asv_name$name_md5 <- str_replace_all(asv_name$name_md5, "_", " ")
  
  # return a dataframe contains taxa ids, full taxonomy, and modified name that can be used in figures
  return(asv_name)
}



####============  Summarize ASV sharing by ASV    ===================####
fun_sum_sharing_by_asv <- function(asv_name, id.list, share_res_asv_v, share_res_asv_f){
  # Input are a dataframe containing ASV taxonomy, and another dataframe containing ASV sharing results for infant*taxa
  
  share_res_asv_v <- filter(share_res_asv_v, row.names(share_res_asv_v) %in% id.list)
  share_res_asv_f <- filter(share_res_asv_f, row.names(share_res_asv_f) %in% id.list)
  
  share_res_byasv1 <- asv_name
  share_res_byasv1$shared = colSums(share_res_asv_v==1)
  share_res_byasv1$not_shared = colSums(share_res_asv_v==4)
  share_res_byasv1$group <- "Maternal vaginal - Infant fecal"
  
  share_res_byasv2 <- asv_name
  share_res_byasv2$shared = colSums(share_res_asv_f==1)
  share_res_byasv2$not_shared = colSums(share_res_asv_f==4)
  share_res_byasv2$group <- "Maternal fecal - Infant fecal"
  
  share_res_byasv <- rbind(share_res_byasv1, share_res_byasv2)
  
  share_res_byasv <- gather(share_res_byasv, type, n_pair, c(shared:not_shared), factor_key=TRUE)
  share_res_byasv$type <- factor(share_res_byasv$type, levels=c("shared", "not_shared"), labels=c("Shared", "Not shared"))
  
  return(share_res_byasv)
  # output is a wide data set: 1 row for each ASV
}




####============  Summarize ASV sharing by ASV at genus or phylum level    ===================####
fun_sum_sharing_by_asv_genusphylum <- function(asv_name, id.list, asv.list, share_res_asv_v, share_res_asv_f){
  
  share_res_asv_v <- filter(share_res_asv_v, row.names(share_res_asv_v) %in% id.list)
  share_res_asv_f <- filter(share_res_asv_f, row.names(share_res_asv_f) %in% id.list)
  
  
  junk1 <- share_res_asv_f
  colnames(junk1) <- sub("^X", "", colnames(junk1))
  junk1 <- as.data.frame(t(junk1))
  junk1 <- merge(asv_name, junk1, by.x=c("asv"), by.y='row.names')
  
  junk1 <- filter(junk1, asv %in% asv.list) # so summarize at genus or phylum level
  
  res_1 <- data.frame(matrix(ncol = 0, nrow = 1))
  res_1$n_shared_any <- sum(apply(junk1, 2, function(col) any(col == 1))) 
  res_1$n_not_shared <- sum(apply(junk1, 2, function(col) all(col != 1) & any(col == 4)))
  res_1$group <- "Maternal fecal - Infant fecal"
  
  
  
  junk2 <- share_res_asv_v
  colnames(junk2) <- sub("^X", "", colnames(junk2))
  junk2 <- as.data.frame(t(junk2))
  junk2 <- merge(asv_name, junk2, by.x=c("asv"), by.y='row.names')
  
  junk2 <- filter(junk2, asv %in% asv.list) # so summarize at genus or phylum level
  
  res_2 <- data.frame(matrix(ncol = 0, nrow = 1))
  res_2$n_shared_any <- sum(apply(junk2, 2, function(col) any(col == 1))) 
  res_2$n_not_shared <- sum(apply(junk2, 2, function(col) all(col != 1) & any(col == 4)))
  res_2$group <- "Maternal vaginal - Infant fecal"
  
  res <- rbind(res_1, res_2)
}


####==============  Conduct Wilcoxon tests  ==============####
fun_wilcoxon_paired_MARCH <- function(dta){
  # convert wide to long
  wide_dta <- spread(dta, type, prop)
  
  # sum maternal contribution
  wide_dta$Mother <- 100 - wide_dta$Other
  
  # test among overall group
  overall_1 <- wilcox.test(wide_dta$Mother, wide_dta$Other, paired=TRUE) # mother vs other
  overall_2 <- wilcox.test(wide_dta$'Maternal fecal', wide_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  # test within stratum of antibiotics
  no_antib_dta <- wide_dta %>% filter(med_antib_pn=="No antibiotics")
  no_antib_1 <- wilcox.test(no_antib_dta$Mother, no_antib_dta$Other, paired=TRUE) # mother vs other
  no_antib_2 <- wilcox.test(no_antib_dta$'Maternal fecal', no_antib_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  antib_dta <- wide_dta %>% filter(med_antib_pn=="Took antibiotics")
  antib_1 <- wilcox.test(antib_dta$Mother, antib_dta$Other, paired=TRUE) # mother vs other
  antib_2 <- wilcox.test(antib_dta$'Maternal fecal', antib_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  # test within stratum of delivery model
  vaginal_dta <- wide_dta %>% filter(birth_delivm=="Vaginal delivery")
  vaginal_1 <- wilcox.test(vaginal_dta$Mother, vaginal_dta$Other, paired=TRUE) # mother vs other
  vaginal_2 <- wilcox.test(vaginal_dta$'Maternal fecal', vaginal_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  no_vaginal_dta <- wide_dta %>% filter(birth_delivm=="Cesarean delivery")
  no_vaginal_1 <- wilcox.test(no_vaginal_dta$Mother, no_vaginal_dta$Other, paired=TRUE) # mother vs other
  no_vaginal_2 <- wilcox.test(no_vaginal_dta$'Maternal fecal', no_vaginal_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  # test within stratum of feeding type
  ebm_dta <- wide_dta %>% filter(feeding=="Exclusive breastmilk")
  ebm_1 <- wilcox.test(ebm_dta$Mother, ebm_dta$Other, paired=TRUE) # mother vs other
  ebm_2 <- wilcox.test(ebm_dta$'Maternal fecal', ebm_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  ef_dta <- wide_dta %>% filter(feeding=="Exclusive formula")
  ef_1 <- wilcox.test(ef_dta$Mother, ef_dta$Other, paired=TRUE) # mother vs other
  ef_2 <- wilcox.test(ef_dta$'Maternal fecal', ef_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  mix_dta <- wide_dta %>% filter(feeding=="Mixed")
  mix_1 <- wilcox.test(mix_dta$Mother, mix_dta$Other, paired=TRUE) # mother vs other
  mix_2 <- wilcox.test(mix_dta$'Maternal fecal', mix_dta$'Maternal vaginal', paired=TRUE) # fecal vs vaginal
  
  
  # summarize p values
  output <- data.frame(group=c("overall", "no antibio", "antibio", "vaginal", "c-section", "ebm", "ef", "mixed"),
                       mo_vs_oth=c(overall_1$p.value, no_antib_1$p.value, antib_1$p.value, vaginal_1$p.value, no_vaginal_1$p.value, ebm_1$p.value, ef_1$p.value, mix_1$p.value),
                       fecal_vs_vaginal=c(overall_2$p.value, no_antib_2$p.value, antib_2$p.value, vaginal_2$p.value, no_vaginal_2$p.value, ebm_2$p.value, ef_2$p.value, mix_2$p.value))
  # keep 4 digits
  output$mo_vs_oth <- round(output$mo_vs_oth, digits=4)
  output$fecal_vs_vaginal <- round(output$fecal_vs_vaginal, digits=4)
  
  
  return(output)
}


fun_wilcoxon_paired_MARCH_2 <- function(wide_dta){

  # test among overall group
  overall_1 <- wilcox.test(wide_dta$sum_not_shared, wide_dta$sum_shared, paired=TRUE) # mother vs other
  overall_2 <- wilcox.test(wide_dta$stool, wide_dta$vagina, paired=TRUE) # fecal vs vaginal
  
  # test within stratum of antibiotics
  no_antib_dta <- wide_dta %>% filter(med_antib_pn=="No antibiotics")
  no_antib_1 <- wilcox.test(no_antib_dta$sum_not_shared, no_antib_dta$sum_shared, paired=TRUE) # mother vs other
  no_antib_2 <- wilcox.test(no_antib_dta$stool, no_antib_dta$vagina, paired=TRUE) # fecal vs vaginal
  
  antib_dta <- wide_dta %>% filter(med_antib_pn=="Took antibiotics")
  antib_1 <- wilcox.test(antib_dta$sum_not_shared, antib_dta$sum_shared, paired=TRUE) # mother vs other
  antib_2 <- wilcox.test(antib_dta$stool, antib_dta$vagina, paired=TRUE) # fecal vs vaginal
  
  # test within stratum of delivery model
  vaginal_dta <- wide_dta %>% filter(birth_delivm=="Vaginal delivery")
  vaginal_1 <- wilcox.test(vaginal_dta$sum_not_shared, vaginal_dta$sum_shared, paired=TRUE) # mother vs other
  vaginal_2 <- wilcox.test(vaginal_dta$stool, vaginal_dta$vagina, paired=TRUE) # fecal vs vaginal
  
  no_vaginal_dta <- wide_dta %>% filter(birth_delivm=="Cesarean delivery")
  no_vaginal_1 <- wilcox.test(no_vaginal_dta$sum_not_shared, no_vaginal_dta$sum_shared, paired=TRUE) # mother vs other
  no_vaginal_2 <- wilcox.test(no_vaginal_dta$stool, no_vaginal_dta$vagina, paired=TRUE) # fecal vs vaginal
  
  # test within stratum of feeding type
  bm_dta <- wide_dta %>% filter(bf=="Yes")
  bm_1 <- wilcox.test(bm_dta$sum_not_shared, bm_dta$sum_shared, paired=TRUE) # mother vs other
  bm_2 <- wilcox.test(bm_dta$stool, bm_dta$vagina, paired=TRUE) # fecal vs vaginal
  
  ef_dta <- wide_dta %>% filter(bf=="No")
  ef_1 <- wilcox.test(ef_dta$sum_not_shared, ef_dta$sum_shared, paired=TRUE) # mother vs other
  ef_2 <- wilcox.test(ef_dta$stool, ef_dta$vagina, paired=TRUE) # fecal vs vaginal
  
  
  # summarize p values
  output <- data.frame(group=c("overall", "no antibio", "antibio", "vaginal", "c-section", "bm", "ef"),
                       mo_vs_oth=c(overall_1$p.value, no_antib_1$p.value, antib_1$p.value, vaginal_1$p.value, no_vaginal_1$p.value, bm_1$p.value, ef_1$p.value),
                       fecal_vs_vaginal=c(overall_2$p.value, no_antib_2$p.value, antib_2$p.value, vaginal_2$p.value, no_vaginal_2$p.value, bm_2$p.value, ef_2$p.value))
  # keep 4 digits
  output$mo_vs_oth <- round(output$mo_vs_oth, digits=4)
  output$fecal_vs_vaginal <- round(output$fecal_vs_vaginal, digits=4)
  
  
  return(output)
}




fun_wilcoxon_paired <- function(dta){
  # convert wide to long
  wide_dta <- spread(dta, type, prop)
  
  # sum maternal contribution
  wide_dta$Mother <- 100 - wide_dta$Other
  
  # test among overall group
  overall_1 <- wilcox.test(wide_dta$Mother, wide_dta$Other, paired=TRUE) # mother vs other

  # test within stratum of antibiotics
  no_antib_dta <- wide_dta %>% filter(med_antib_pn=="No antibiotics")
  no_antib_1 <- wilcox.test(no_antib_dta$Mother, no_antib_dta$Other, paired=TRUE) # mother vs other

  antib_dta <- wide_dta %>% filter(med_antib_pn=="Took antibiotics")
  antib_1 <- wilcox.test(antib_dta$Mother, antib_dta$Other, paired=TRUE) # mother vs other

  # test within stratum of delivery model
  vaginal_dta <- wide_dta %>% filter(birth_delivm=="Vaginal delivery")
  vaginal_1 <- wilcox.test(vaginal_dta$Mother, vaginal_dta$Other, paired=TRUE) # mother vs other

  no_vaginal_dta <- wide_dta %>% filter(birth_delivm=="Cesarean delivery")
  no_vaginal_1 <- wilcox.test(no_vaginal_dta$Mother, no_vaginal_dta$Other, paired=TRUE) # mother vs other

  # summarize p values
  output <- data.frame(group=c("overall", "no antibio", "antibio", "vaginal", "c-section"),
                       mo_vs_oth=c(overall_1$p.value, no_antib_1$p.value, antib_1$p.value, vaginal_1$p.value, no_vaginal_1$p.value))
  
  # keep 4 digits
  output$mo_vs_oth <- round(output$mo_vs_oth, digits=4)
  
  return(output)
}


####==============  Run mixed effect models for all cohorts together  ==================####
# Allow for both random intercepts and random slopes
run_lme_by_cohort <- function(df, x, y, group, cohort, correlation_structure = NULL) {
  
  # Initialize lists to store the results
  combined_outs <- list()
  combined_dfs <- list()
  
  # Get the unique cohort names
  cohorts <- unique(df[[cohort]])
  
  # Loop through each cohort
  for(i in seq_along(cohorts)) {
    
    # Subset the data for this cohort
    df_cohort <- df[df[[cohort]] == cohorts[i], ]
    
    # Fit the lme model
    model <- lme(contribution ~ age_mo, random = ~ 1 + age_mo | ParticipantID,
                 correlation = correlation_structure, # define correlation structure
                 control = lmeControl(opt = "optim", optCtrl = list(maxit = 10000)),
                 data = df_cohort)
    
    # Store the model summary for the fixed effects
    model_out <- as.data.frame(t(summary(model)$tTable[2,]))
    model_out$AIC <- AIC(model)
    model_out$BIC <- BIC(model)
    
    model_out[[cohort]] <- cohorts[i]
    combined_out <- model_out
    combined_outs[[cohorts[i]]] <- combined_out
    
    
    # Create a new data frame for prediction
    newdata <- data.frame(x = seq(min(df_cohort[[x]]), max(df_cohort[[x]])))
    names(newdata) <- x
    
    # Predict y values
    newdata$y_pred <- predict(model, level=0, newdata)
    
    newdata[[cohort]] <- cohorts[i]
    combined_df <- newdata
    
    # Store the combined data
    combined_dfs[[cohorts[i]]] <- combined_df
  }
  
  return(list(combined_outs = combined_outs, combined_dfs = combined_dfs))
}


####==============  Create Venn diagram to show overlap between infant, maternal fecal, and maternal vaginal microbiota  ==================####
# Note, this code creates non-proportional venn diagrams
venn_3group <- function(share_res_byinfant_wide){
  # Summarize mean sharing count
  both_avg = mean(share_res_byinfant_wide$both)
  vagina_avg = mean(share_res_byinfant_wide$vagina)
  stool_avg = mean(share_res_byinfant_wide$stool)
  
  infant_only_avg = mean(share_res_byinfant_wide$sum_infant - share_res_byinfant_wide$both - share_res_byinfant_wide$vagina - share_res_byinfant_wide$stool)
  vagina_only_avg = mean(share_res_byinfant_wide$sum_vagina - share_res_byinfant_wide$overlap - share_res_byinfant_wide$vagina)
  stool_only_avg = mean(share_res_byinfant_wide$sum_stool - share_res_byinfant_wide$overlap - share_res_byinfant_wide$stool)
  
  
  overlap_adj_avg = mean(share_res_byinfant_wide$overlap - share_res_byinfant_wide$both)
  
  # Define a sequence for each of the categories and overlaps
  infant_seq = 1:infant_only_avg
  stool_seq = (max(infant_seq)+1):(max(infant_seq)+stool_only_avg)
  vaginal_seq = (max(stool_seq)+1):(max(stool_seq)+vagina_only_avg)
  
  both_seq = (max(vaginal_seq)+1):(max(vaginal_seq)+both_avg)
  vagina_seq = (max(both_seq)+1):(max(both_seq)+vagina_avg)
  stool_seq_overlap = (max(vagina_seq)+1):(max(vagina_seq)+stool_avg)
  fecal_vaginal_seq = (max(stool_seq_overlap)+1):(max(stool_seq_overlap)+overlap_adj_avg)
  
  # Define the sets
  sets = list(
    Infant = c(infant_seq, both_seq, vagina_seq, stool_seq_overlap),
    MaternalFecal = c(stool_seq, both_seq, stool_seq_overlap, fecal_vaginal_seq),
    MaternalVaginal = c(vaginal_seq, both_seq, vagina_seq, fecal_vaginal_seq)
  )
  
  # Plot the Venn diagram
  venn <- Venn(sets)
  data <- process_data(venn)
  
  n <- sum(venn_region(data)$count)
  
  p <- ggplot() +
    geom_sf(aes(fill=factor(name, 
                               levels=c("Infant..MaternalVaginal", "MaternalFecal..MaternalVaginal",
                                        "Infant..MaternalFecal..MaternalVaginal", "Infant..MaternalFecal", 
                                        "MaternalVaginal", "Infant", "MaternalFecal"))), data = venn_region(data)) +
    geom_sf(lwd=0.8, color = "black", data = venn_setedge(data), show.legend = F) +
    geom_text(aes(x=c(250, 500, 750), y=c(100, 910, 100), label = factor(name, levels=c("Infant", "MaternalFecal", "MaternalVaginal"),
                                                                         labels=c("Infant fecal", "Maternal fecal", "Maternal vaginal"))), 
              fontface="bold", size=7, data = venn_setlabel(data)) +
    geom_sf_text(aes(label=paste0(count, "\n(", round(count/n*100, digits = 1), "%)")), size=5.5, data = venn_region(data)) +
    scale_x_continuous(expand = expansion(mult = .2)) + 
    theme_void() + 
    scale_fill_manual(values=c("#E69F00", "#CC79A7", "#90C987", "#0072B2", "white","white", "white")) +
    theme(legend.position = "none") 
}


####===============   Compare number of shared ASVs from each maternal source between groups of variables   ===============####
test_ASV_byvariabels <- function(wide_data){
  
  wide_data$sum_share_stool <- wide_data$stool + wide_data$both
  wide_data$sum_share_vagina <- wide_data$vagina + wide_data$both
  
  # convert from wide to long
  data <- gather(wide_data, type, n_count, c(sum_share_stool:sum_share_vagina, sum_not_shared), factor_key=TRUE)
  
  data$type <- factor(data$type, 
                      levels=c("sum_not_shared","sum_share_stool","sum_share_vagina"),
                      labels=c("Not shared with either maternal vaginal or fecal microbiota", 
                               "Shared with maternal fecal microbiota", "Shared with maternal vaginal microbiota"))
  # define the color I want
  color <- c("#BBBBBB","#00507D", "#A16F00")
  
  
  # by prenatal antibiotics
  stat.test <- data %>%
    group_by(type) %>%
    wilcox_test(n_count ~ med_antib_pn) %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "med_antib_pn")
  stat.test$xmax <- stat.test$xmax - 1
  stat.test$xmin <- stat.test$xmin - 1
  
  p1 <- 
    data %>%
    ggplot(aes(y=n_count, x=factor(med_antib_pn, levels=c("No antibiotics", "Took antibiotics"), labels=c("No", "Yes")))) +
    geom_boxplot(lwd=0.8, width=0.8, aes(color=type), outlier.shape = NA) + 
    geom_point(size=3, aes(color=type), position=position_jitter(w = 0.1, h = 0.1), alpha=0.5) +
    ylab("Number of ASV")  + xlab("") + labs(color="Sharing", title="A. Prenatal Antibiotics") + 
    facet_grid(~type) +
    theme_bw() +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=15, face="bold"), axis.ticks.x = element_blank(),
          axis.title.x=element_blank(), axis.title = element_text(size=15, face="bold"),
          legend.title = element_text(size=15, face="bold"), legend.text = element_text(size=15), legend.direction = "horizontal",
          legend.position = "bottom",
          panel.grid = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(),
          panel.border = element_blank(), 
          axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour="black"),
          plot.title = element_text(size=17, face="bold")) +
    scale_color_manual(values=color) + 
    stat_summary(fun.y="mean", aes(color=type), shape=18, size=1.7) +
    stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size=5) +
    scale_y_continuous(limits=c(0, 280), breaks=c(0, 50, 100, 150, 200, 250))
  
  
  # by delivery mode
  stat.test <- data %>%
    group_by(type) %>%
    wilcox_test(n_count ~ birth_delivm) %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "birth_delivm")
  stat.test$xmax <- stat.test$xmax - 1
  stat.test$xmin <- stat.test$xmin - 1
  
  p2 <- 
    data %>%
    ggplot(aes(y=n_count, x=factor(birth_delivm, levels=c("Vaginal delivery", "Cesarean delivery"), labels=c("Vaginal", "C-section")))) +
    geom_boxplot(lwd=0.8, width=0.8, aes(color=type), outlier.shape = NA) + 
    geom_point(size=3, aes(color=type), position=position_jitter(w = 0.1, h = 0.1), alpha=0.5) +
    ylab("Number of ASV")  + xlab("") + labs(color="Sharing", title="B. Delivery Mode") + 
    facet_grid(~type) +
    theme_bw() +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=15, face="bold"), axis.ticks.x = element_blank(),
          axis.title.x=element_blank(), axis.title = element_text(size=15, face="bold"),
          legend.title = element_text(size=15, face="bold"), legend.text = element_text(size=15), legend.direction = "horizontal",
          legend.position = "bottom",
          panel.grid = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(),
          panel.border = element_blank(), 
          axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour="black"),
          plot.title = element_text(size=17, face="bold")) +
    scale_color_manual(values=color) + 
    stat_summary(fun.y="mean", aes(color=type), shape=18, size=1.7) +
    stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size=5) +
    scale_y_continuous(limits=c(0, 280), breaks=c(0, 50, 100, 150, 200, 250))
  p2
  
  # by breastfeeding
  stat.test <- data[data$bf!="Missing",] %>%
    group_by(type) %>%
    wilcox_test(n_count ~ bf) %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "bf")
  stat.test$xmax <- stat.test$xmax - 1
  stat.test$xmin <- stat.test$xmin - 1
  
  p3 <- 
    data[data$bf!="Missing",]  %>%
    ggplot(aes(y=n_count, x=factor(bf, levels=c("Yes", "No")))) +
    geom_boxplot(lwd=0.8, width=0.8, aes(color=type), outlier.shape = NA) + 
    geom_point(size=3, aes(color=type), position=position_jitter(w = 0.1, h = 0.1), alpha=0.5) +
    ylab("Number of ASV")  + xlab("") + labs(color="Sharing", title="C. Ever Breastfed") + 
    facet_grid(~type) +
    theme_bw() +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=15, face="bold"), axis.ticks.x = element_blank(),
          axis.title.x=element_blank(), axis.title = element_text(size=15, face="bold"),
          legend.title = element_text(size=15, face="bold"), legend.text = element_text(size=15), legend.direction = "horizontal",
          legend.position = "bottom",
          panel.grid = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(),
          panel.border = element_blank(), 
          axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour="black"),
          plot.title = element_text(size=17, face="bold")) +
    scale_color_manual(values=color) + 
    stat_summary(fun.y="mean", aes(color=type), shape=18, size=1.7) +
    stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size=5) +
    scale_y_continuous(limits=c(0, 280), breaks=c(0, 50, 100, 150, 200, 250))
  
  p <- 
  ggarrange(p1 + ylab("Number of ASV") + theme(legend.direction = "horizontal"), 
            p2 + theme(axis.title.y = element_blank()), 
            p3 + theme(axis.title.y = element_blank()),
            widths = c(1, 0.95, 0.95), nrow=1, common.legend = TRUE, legend="bottom")
  return(p)
}


####===============   Compare FEAST estimated contributions from each maternal source between groups of variables   ===============####
test_FEAST_byvariabels <- function(march){
  # define the color I want
  color <- c("#BBBBBB","#00507D", "#A16F00")
  
  
  stat.test <- march %>%
    group_by(type) %>%
    wilcox_test(prop ~ med_antib_pn) %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "med_antib_pn")
  
  p1 <-
    march  %>%
    ggplot(aes(y=prop, x=factor(med_antib_pn, levels=c("No", "Yes")))) +
    geom_boxplot(lwd=0.8, width=0.8, aes(color=type), outlier.shape = NA) + geom_point(size=3, aes(color=type), position=position_jitter(w = 0.1, h = 0.1), alpha=0.5) +
    ylab("Contribution of Source (%)")  + xlab("") + labs(color="Source") + ggtitle("A. Prenatal Antibiotic Use") + 
    facet_grid(~type) +
    theme_bw() +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=15, face="bold"), axis.ticks.x = element_blank(),
          axis.title.x=element_blank(), axis.title = element_text(size=15, face="bold"),
          legend.title = element_text(size=15, face="bold"), legend.text = element_text(size=15), legend.direction = "vertical",
          panel.grid = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(),
          panel.border = element_blank(), 
          axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour="black"),
          plot.title = element_text(size=17, face="bold")) +
    scale_color_manual(values=color) + stat_summary(fun.y="mean", aes(color=type), shape=18, size=1.7) +
    stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size=5) +
    scale_y_continuous(limits=c(0, 110), breaks=c(0, 25, 50, 75, 100))
  p1
  
  
  stat.test <- march %>%
    group_by(type) %>%
    wilcox_test(prop ~ birth_delivm) %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "birth_delivm")
  
  p2 <-
    march  %>%
    ggplot(aes(y=prop, x=factor(birth_delivm, levels=c("Vaginal", "C-section")))) +
    geom_boxplot(lwd=0.8, width=0.8, aes(color=type), outlier.shape = NA) + geom_point(size=3, aes(color=type), position=position_jitter(w = 0.1, h = 0.1), alpha=0.5) +
    ylab("Contribution of Source (%)")  + xlab("") + labs(color="Source", title="B. Delivery Mode") + 
    facet_grid(~type) +
    theme_bw() +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=15, face="bold"), axis.ticks.x = element_blank(),
          axis.title.x=element_blank(), axis.title = element_text(size=15, face="bold"),
          legend.title = element_text(size=15, face="bold"), legend.text = element_text(size=15), legend.direction = "vertical",
          panel.grid = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(),
          panel.border = element_blank(), 
          axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour="black"),
          plot.title = element_text(size=17, face="bold")) +
    scale_color_manual(values=color) + stat_summary(fun.y="mean", aes(color=type), shape=18, size=1.7) +
    stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size=5) +
    scale_y_continuous(limits=c(0, 110), breaks=c(0, 25, 50, 75, 100))
  p2
  
  
  stat.test <- march %>%
    group_by(type) %>%
    wilcox_test(prop ~ bf) %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "bf")

  p3 <- march  %>%
    ggplot(aes(y=prop, x=bf)) +
    geom_boxplot(lwd=0.8, width=0.8, aes(color=type), outlier.shape = NA) + geom_point(size=3, aes(color=type), position=position_jitter(w = 0.1, h = 0.1), alpha=0.5) +
    ylab("Contribution of Source (%)")  + xlab("") + labs(color="Source", title="C. Ever Breastfed") + 
    facet_grid(~type) +
    theme_bw() +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=15, face="bold"), axis.ticks.x = element_blank(),
          axis.title.x=element_blank(), axis.title = element_text(size=15, face="bold"),
          legend.title = element_text(size=15, face="bold"), legend.text = element_text(size=15), legend.direction = "vertical",
          panel.grid = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(),
          panel.border = element_blank(), 
          axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour="black"),
          plot.title = element_text(size=17, face="bold")) +
    scale_color_manual(values=color) + stat_summary(fun.y="mean", aes(color=type), shape=18, size=1.7) +
    stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size=5) +
    scale_y_continuous(limits=c(0, 110), breaks=c(0, 25, 50, 75, 100))
  
  p3
  
  p <- ggarrange(p1 + ylab("Contribution of Source (%)") + theme(legend.direction = "horizontal"), 
                 p2 + theme(axis.title.y = element_blank()), 
                 p3 + theme(axis.title.y = element_blank()),
                 widths = c(1, 0.95, 0.95), nrow=1, common.legend = TRUE, legend="bottom")
  return(p)
  
}



####=============   Plot Bifidobacterium and Bacteroides and compare sharing frequencies between groups of variables    ============####
Bifido_Bac_plot <- function(asv_dta){
# input is either Bifidobacterium or Bacteroides data with clean asv name 
  
  # define color that will be used
  my_color <- c("white","#D55E00") # white for not shared, red for shared

  # only keep ASVs with evidence of sharing
  junk <- asv_dta[asv_dta$type=="Shared" & asv_dta$variable=="Overall",] %>%
    group_by(name) %>%
    mutate(sum_any_share=sum(n_pair))
  junk <- junk[junk$sum_any_share>0, ]
  
  # calculate pairs of any share overall, which will be used to sort ASVs
  asv_dta <- filter(asv_dta, asv %in% as.list(junk[junk$sum_any_share>0, ]$asv))
  asv_dta <- merge(asv_dta, junk[, c("asv", "sum_any_share", "group")], by=c("asv", "group"))
  
  
  # calculate proportion of shared, which will be added as labels in the figure
  asv_dta <- asv_dta %>%
    group_by(name, group, variable) %>%
    mutate(sum_n=sum(n_pair), label=ifelse(type=="Not shared" | sum_n==0 | n_pair==0, "", paste0(round(n_pair/sum_n*100, digits = 1), "%")))
  
  
  ## fig 0: overall
  fig_dta <- asv_dta[asv_dta$variable=="Overall",]
  
  fig_overall <-
    ggplot(fig_dta, aes(x=variable, y=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black") +  
    geom_text(data=fig_dta, aes(x=variable, y=n_pair, label=label, group=type), size = 3.5, position = position_stack(vjust=0.5)) +  
    facet_grid(reorder(name_md5, sum_any_share, decreasing=TRUE)~group, scales = "free_y", switch="y", space = "free_y") +
    ylab("Number of Mother-Infant Dyads") + xlab("Overall") + labs(fill="") +
    theme_bw() +
    theme(axis.title.x = element_text(size=15, face="bold"), axis.title.y = element_text(size=14), axis.ticks.y=element_blank(),
          axis.text.x = element_text(size=13), axis.text.y = element_blank(),
          legend.position="bottom", legend.direction="horizontal", legend.text = element_text(size=13),
          legend.title = element_text(size=13, face="bold"),
          panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust=1), strip.placement = "outside", 
          strip.background.y = element_blank(), strip.text.y = element_text(size=15, face="italic"), 
          strip.text.x = element_text(size=15, color="white"), strip.background.x = element_rect(color=NA, fill="#999999"),
          panel.border = element_blank() , axis.line.x = element_line(colour = "black")) +
    coord_flip() + ylim(0, 25) +
    scale_fill_manual(values=my_color) + scale_x_discrete(position="top") 
  
  
  
  ## fig 1: by prenatal antibiotics
  fig_dta <- asv_dta[asv_dta$variable=="No antibiotics" | asv_dta$variable=="Took antibiotics",]
  fig_dta$variable <- factor(fig_dta$variable, levels=c("Took antibiotics", "No antibiotics"), labels=c("Yes", "No"))
  
    # conduct Fisher's test
  test_data <- fig_dta[, 1:17] %>%
    group_by(group) %>%
    spread(key = type, value = n_pair)
  
  test_result <- test_data %>%
    group_by(name_md5, group) %>%
    do(test = fisher.test(matrix(c(.$Shared[.$variable == "No"], 
                                   .$Shared[.$variable == "Yes"], 
                                   .$'Not shared'[.$variable == "No"], 
                                   .$'Not shared'[.$variable == "Yes"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  
    # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  fig_dta <- merge(fig_dta, test_result[, c("name_md5", "group",  "significance")], by = c("name_md5", "group"), all.x = TRUE)
  
  # make plot
  fig_antib <-
    ggplot(fig_dta, aes(x=variable, y=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black") +  
    geom_text(data=fig_dta, aes(x=variable, y=n_pair, label=label, group=type), size = 3.5, position = position_stack(vjust=0.5)) +  
    facet_grid(reorder(name_md5, sum_any_share, decreasing=TRUE)~group, scales = "free_y", switch="y", space = "free_y") +
    ylab("Number of Mother-Infant Dyads") + xlab("Prenatal Antibiotic Use") + labs(fill="") +
    theme_bw() +
    theme(axis.title.x = element_text(size=15, face="bold"), axis.title.y = element_text(size=14), axis.ticks.y=element_blank(),
          axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
          legend.position="bottom", legend.direction="horizontal", legend.text = element_text(size=13),
          legend.title = element_text(size=13, face="bold"),
          panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust=1), strip.placement = "outside", 
          strip.background.y = element_blank(), strip.text.y = element_text(size=15, face="italic"), 
          strip.text.x = element_text(size=15, color="white"), strip.background.x = element_rect(color=NA, fill="#999999"),
          panel.border = element_blank() , axis.line.x = element_line(colour = "black")) +
    coord_flip() + ylim(0, 13) +
    scale_fill_manual(values=my_color) + scale_x_discrete(position="top") +
    geom_text(aes(x = 1, y = 13, label = significance), size=9)
  
  
  ## fig 2: by delivery mode
  fig_dta <- asv_dta[asv_dta$variable=="C-section" | asv_dta$variable=="Vaginal",]
  fig_dta$variable <- factor(fig_dta$variable, levels=c("C-section", "Vaginal"))
  
  # conduct Fisher's test
  test_data <- fig_dta[, 1:17] %>%
    group_by(group) %>%
    spread(key = type, value = n_pair)
  
  test_result <- test_data %>%
    group_by(name_md5, group) %>%
    do(test = fisher.test(matrix(c(.$Shared[.$variable == "Vaginal"], 
                                   .$Shared[.$variable == "C-section"], 
                                   .$'Not shared'[.$variable == "Vaginal"], 
                                   .$'Not shared'[.$variable == "C-section"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  fig_dta <- merge(fig_dta, test_result[, c("name_md5", "group",  "significance")], by = c("name_md5", "group"), all.x = TRUE)
  
  # make plot
  fig_delivery <-
    ggplot(fig_dta, aes(x=variable, y=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black") +  
    geom_text(data=fig_dta, aes(x=variable, y=n_pair, label=label, group=type), size = 3.5, position = position_stack(vjust=0.5)) +  
    facet_grid(reorder(name_md5, sum_any_share, decreasing=TRUE)~group, scales = "free_y", switch="y", space = "free_y") +
    ylab("Number of Mother-Infant Dyads") + xlab("Delivery Mode") + labs(fill="") +
    theme_bw() +
    theme(axis.title.x = element_text(size=15, face="bold"), axis.title.y = element_text(size=14), axis.ticks.y=element_blank(),
          axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
          legend.position="bottom", legend.direction="horizontal", legend.text = element_text(size=13),
          legend.title = element_text(size=13, face="bold"),
          panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust=1), strip.placement = "outside", 
          strip.background.y = element_blank(), strip.text.y = element_text(size=15, face="italic"), 
          strip.text.x = element_text(size=15, color="white"), strip.background.x = element_rect(color=NA, fill="#999999"),
          panel.border = element_blank() , axis.line.x = element_line(colour = "black")) +
    coord_flip() + ylim(0, 20) +
    scale_fill_manual(values=my_color) + scale_x_discrete(position="top") +
    geom_text(aes(x = 1, y = 20, label = significance), size=9)
  
  
  
  ## fig 3: by ever breastfed
  fig_dta <- asv_dta[asv_dta$variable=="Yes" | asv_dta$variable=="No",]
  
    # conduct Fisher's test
  test_data <- fig_dta[, 1:17] %>%
    group_by(group) %>%
    spread(key = type, value = n_pair)
  
  test_result <- test_data %>%
    group_by(name_md5, group) %>%
    do(test = fisher.test(matrix(c(.$Shared[.$variable == "Yes"], 
                                   .$Shared[.$variable == "No"], 
                                   .$'Not shared'[.$variable == "Yes"], 
                                   .$'Not shared'[.$variable == "No"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  
    # Convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  fig_dta <- merge(fig_dta, test_result[, c("name_md5", "group",  "significance")], by = c("name_md5", "group"), all.x = TRUE)
  
    # make plot
  fig_bf <-
    ggplot(fig_dta, aes(x=variable, y=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black") +  
    geom_text(data=fig_dta, aes(x=variable, y=n_pair, label=label, group=type), size = 3.5, position = position_stack(vjust=0.5)) +  
    facet_grid(reorder(name_md5, sum_any_share, decreasing=TRUE)~group, scales = "free_y", switch="y", space = "free_y") +
    ylab("Number of Mother-Infant Dyads") + xlab("Ever Breastfed") + labs(fill="") +
    theme_bw() +
    theme(axis.title.x = element_text(size=15, face="bold"), axis.title.y = element_text(size=14), axis.ticks.y=element_blank(),
          axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
          legend.position="bottom", legend.direction="horizontal", legend.text = element_text(size=13),
          legend.title = element_text(size=13, face="bold"),
          panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust=1), strip.placement = "outside", 
          strip.background.y = element_blank(), strip.text.y = element_text(size=15, face="italic"), 
          strip.text.x = element_text(size=15, color="white"), strip.background.x = element_rect(color=NA, fill="#999999"),
          panel.border = element_blank() , axis.line.x = element_line(colour = "black")) +
    coord_flip() +
    scale_fill_manual(values=my_color) + scale_x_discrete(position="top") +
    geom_text(aes(x = 1, y = 20, label = significance), size=9)
  
  
  # combine plot
  p <-
    ggarrange(fig_overall, fig_antib + theme(strip.text.y=element_text(size=0, color="transparent")), 
            fig_delivery + theme(strip.text.y=element_text(size=0, color="transparent")), 
            fig_bf + theme(strip.text.y=element_text(size=0, color="transparent")), 
            nrow=1, common.legend = TRUE, legend = "bottom", widths = c(0.95, 0.75, 0.8, 0.75))
  
  return(p)
}


####=============   Plot Bifidobacterium and Bacteroides and compare sharing frequencies between groups of variables: Heatmap    ============####
Bifido_Bac_heat <- function(asv_dta){
# input is either Bifidobacterium or Bacteroides data with clean asv name 

  # only keep ASVs with evidence of sharing
  junk <- asv_dta[asv_dta$type=="Shared" & asv_dta$variable=="Overall",] %>%
    group_by(name) %>%
    mutate(sum_any_share=sum(n_pair))
  junk <- junk[junk$sum_any_share>0, ]
  
  # calculate pairs of any share overall, which will be used to sort ASVs
  asv_dta <- filter(asv_dta, asv %in% as.list(junk[junk$sum_any_share>0, ]$asv))
  asv_dta <- merge(asv_dta, junk[, c("asv", "sum_any_share", "group")], by=c("asv", "group"))
  
  
  # calculate proportion of shared, which will be added as labels in the figure
  asv_dta <- asv_dta %>%
    group_by(name, group, variable) %>%
    mutate(sum_n=sum(n_pair), prop=n_pair/sum_n*100,
           label=ifelse(sum_n==0, "NA", paste0(round(n_pair/sum_n*100, digits = 1))))
  
  
  ## fig 0: overall
  fig_dta <- asv_dta[asv_dta$variable=="Overall" & asv_dta$type=="Shared",]
  
  fig_overall_1 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal fecal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile() + coord_fixed() + ylab("Overall") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title = element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), 
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), 
          strip.background.y = element_blank(), strip.text.y = element_text(size=8, face="italic"),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#0072B2") + scale_y_discrete(position = "right")

  fig_overall_2 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal vaginal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile() + coord_fixed() + ylab("Overall") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title.x = element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), 
          axis.title.y = element_text(size=8),
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), 
          strip.background.y = element_blank(), strip.text.y = element_text(size=0, face="italic", color="transparent"),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#E69F00") + scale_y_discrete(position = "right") 

  
  ## fig 1: by prenatal antibiotics
  fig_dta <- asv_dta[asv_dta$variable=="No antibiotics" | asv_dta$variable=="Took antibiotics",]
  fig_dta$variable <- factor(fig_dta$variable, levels=c("Took antibiotics", "No antibiotics"), labels=c("Yes", "No"))
  
  # conduct Fisher's test
  test_data <- fig_dta[, 1:17] %>%
    group_by(group) %>%
    spread(key = type, value = n_pair)
  
  test_result <- test_data %>%
    group_by(name_md5, group) %>%
    do(test = fisher.test(matrix(c(.$Shared[.$variable == "No"], 
                                   .$Shared[.$variable == "Yes"], 
                                   .$'Not shared'[.$variable == "No"], 
                                   .$'Not shared'[.$variable == "Yes"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  fig_dta <- merge(fig_dta, test_result[, c("name_md5", "group",  "significance")], by = c("name_md5", "group"), all.x = TRUE)
  
  fig_dta <- fig_dta[fig_dta$type=="Shared",]
  
  fig_antib_1 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal fecal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile(color="white") + coord_fixed() + ylab("Prenatal Antibiotic Use") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title = element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), 
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), 
          strip.background.y = element_blank(), strip.text.y = element_text(size=0, face="italic", color="transparent"),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#0072B2", na.value="#BBBBBB") + scale_y_discrete(position = "right") +
    geom_text(aes(x = 1, y = 1.5, label = significance), size=3, color="red")
  
  fig_antib_2 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal vaginal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile(color="white") + coord_fixed() + ylab("Prenatal Antibiotic Use") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title.x = element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size=8), axis.text.y = element_text(size=6), 
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), 
          strip.background.y = element_blank(), strip.text.y = element_text(size=0, face="italic", color="transparent"),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#E69F00", na.value="#BBBBBB") + scale_y_discrete(position = "right") +
    geom_text(aes(x = 1, y = 1.5, label = significance), size=3, color="red")
  fig_antib_2
  
  
  ## fig 2: by delivery mode
  fig_dta <- asv_dta[(asv_dta$variable=="C-section" | asv_dta$variable=="Vaginal"),]
  fig_dta$variable <- factor(fig_dta$variable, levels=c("C-section", "Vaginal"))
  
  # conduct Fisher's test
  test_data <- fig_dta[, 1:17] %>%
    group_by(group) %>%
    spread(key = type, value = n_pair)
  
  test_result <- test_data %>%
    group_by(name_md5, group) %>%
    do(test = fisher.test(matrix(c(.$Shared[.$variable == "Vaginal"], 
                                   .$Shared[.$variable == "C-section"], 
                                   .$'Not shared'[.$variable == "Vaginal"], 
                                   .$'Not shared'[.$variable == "C-section"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  fig_dta <- merge(fig_dta, test_result[, c("name_md5", "group",  "significance")], by = c("name_md5", "group"), all.x = TRUE)
  
  fig_dta <- fig_dta[fig_dta$type=="Shared",]
  
  fig_delivery_1 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal fecal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile(color="white") + coord_fixed() + ylab("Delivery Mode") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title = element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), 
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), 
          strip.background.y = element_blank(), strip.text.y = element_text(size=0, face="italic", color="transparent"),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#0072B2", na.value="#BBBBBB") + scale_y_discrete(position = "right") +
    geom_text(aes(x = 1, y = 1.5, label = significance), size=3, color="red")
  
  fig_delivery_2 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal vaginal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile(color="white") + coord_fixed() + ylab("Delivery Mode") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title.x = element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size=8), axis.text.y = element_text(size=6), 
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), 
          strip.background.y = element_blank(), strip.text.y = element_text(size=0, face="italic", color="transparent"),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#E69F00", na.value="#BBBBBB") + scale_y_discrete(position = "right") +
    geom_text(aes(x = 1, y = 1.5, label = significance), size=3, color="red")
  
  
  ## fig 3: by ever breastfed
  fig_dta <- asv_dta[asv_dta$variable=="Yes" | asv_dta$variable=="No",]
  
  # conduct Fisher's test
  test_data <- fig_dta[, 1:17] %>%
    group_by(group) %>%
    spread(key = type, value = n_pair)
  
  test_result <- test_data %>%
    group_by(name_md5, group) %>%
    do(test = fisher.test(matrix(c(.$Shared[.$variable == "Yes"], 
                                   .$Shared[.$variable == "No"], 
                                   .$'Not shared'[.$variable == "Yes"], 
                                   .$'Not shared'[.$variable == "No"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  
  # Convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  fig_dta <- merge(fig_dta, test_result[, c("name_md5", "group",  "significance")], by = c("name_md5", "group"), all.x = TRUE)
  
  fig_dta <- fig_dta[fig_dta$type=="Shared",]
  
  fig_bf_1 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal fecal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile(color="white") + coord_fixed() + ylab("Ever Breastfed") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title = element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), 
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), 
          strip.background.y = element_blank(), strip.text.y = element_text(size=0, face="italic", color="transparent"),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#0072B2", na.value="#BBBBBB") + scale_y_discrete(position = "right") +
    geom_text(aes(x = 1, y = 1.5, label = significance), size=3, color="red")
  
  fig_bf_2 <-
    ggplot(fig_dta[fig_dta$group=="Sharing with maternal vaginal microbiota",], aes(x=1, y=variable, fill=prop)) +
    geom_tile(color="white") + coord_fixed() + ylab("Ever Breastfed") + labs(fill="") +
    geom_text(aes(label = label), color = "black", size = 2) +
    facet_wrap(~ reorder(name_md5, sum_any_share, decreasing=TRUE), ncol = 1, strip.position = "left") + 
    theme_bw() +   
    theme(axis.title.x = element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size=8), axis.text.y = element_text(size=6), 
          legend.position="none", 
          panel.grid = element_blank(),
          strip.placement = "outside", strip.text = element_blank(), 
          strip.background = element_blank(),
          panel.border = element_blank()) +
    scale_fill_gradient(limits=c(0, 100), low = "white", high = "#E69F00", na.value="#BBBBBB") + scale_y_discrete(position = "right")  +
    geom_text(aes(x = 1, y = 1.5, label = significance), size=3, color="red")
  
  ## arrange together
  p1 <- ggarrange(fig_overall_1 + theme(plot.margin=unit(c(0, 0, 0, 0),"cm")), 
                  fig_overall_2 + theme(plot.margin=unit(c(0, 0, 0, -2.5),"cm")), nrow=1, widths=c(1.2, 0.7))
  
  p2 <- ggarrange(fig_antib_1 + theme(plot.margin=unit(c(0, 0, 0, 0),"cm")), 
                  fig_antib_2 + theme(plot.margin=unit(c(0, 0, 0, -2.5),"cm")), nrow=1, widths=c(0.3, 0.5))
  
  p3 <- ggarrange(fig_delivery_1 + theme(plot.margin=unit(c(0, 0, 0, 0),"cm")), 
                  fig_delivery_2 + theme(plot.margin=unit(c(0, 0, 0, -2.5),"cm")), nrow=1, widths=c(0.3, 0.7))
  
  p4 <-  ggarrange(fig_bf_1 + theme(plot.margin=unit(c(0, 0, 0, 0),"cm")), 
                   fig_bf_2 + theme(plot.margin=unit(c(0, 0, 0, -2.5),"cm")), nrow=1, widths=c(0.3, 0.5))
  
  ## combine plot
  p <- 
    ggarrange(p1 + theme(plot.margin=unit(c(0, -2.2, 0, 0), "cm")), 
              p2 + theme(plot.margin=unit(c(0, -2, 0, 0), "cm")), 
              p3 + theme(plot.margin=unit(c(0, -2.2, 0, 0), "cm")), 
              p4 + theme(plot.margin=unit(c(0, -2, 0, 0), "cm")), 
              nrow=1, widths=c(1.5, 0.6, 0.7, 0.65))  
  return(p)
}


####=============   Plot how many ASVs of Bifidobacterium and Bacteroides genus are shared, overall and between groups of variables    ============####
Bifido_Bac_genus <- function(asv_dta){
  asv_dta <- asv_dta %>%
    group_by(name, group, variable) %>%
    mutate(sum_n=sum(n_pair))
  
  asv_dta <- asv_dta[asv_dta$type=="Shared",] %>% 
    group_by(group, variable) %>%
    summarise(n_detectable=sum(sum_n>0), n_shared=sum(n_pair>0), n_not_shared=sum(sum_n>0) - sum(n_pair>0))
  
  uplimt <- max(asv_dta$n_detectable) + 1
  
  fig_dta <- gather(asv_dta, type, n_pair, n_shared:n_not_shared, factor_key=TRUE)
  fig_dta$prop <- fig_dta$n_pair / fig_dta$n_detectable * 100
  fig_dta$type <- factor(fig_dta$type, levels=c("n_not_shared", "n_shared"), labels=c("Not shared", "Shared"))
  
  
  # overall
  plot_dta <- filter(fig_dta, variable=="Overall")
  
  p_overall_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota", ] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text = element_blank(),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) 
  
  p_overall_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota", ] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text = element_blank(),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank())
  
  
  p_overall <- ggarrange(p_overall_1 + theme(plot.margin=unit(c(0.5, 0, 0, 0),"cm")), 
                         p_overall_2 + theme(plot.margin=unit(c(-0.2,0, 0.5, 0),"cm")), 
                         nrow=2, ncol=1, heights = c(1, 0.9))
  p_overall <- annotate_figure(p_overall, right=text_grob("Overall",  size = 8, rot = -90))
  
  
  # by prenatal antibiotic use
  test_result <- asv_dta[asv_dta$variable=="No antibiotics" | asv_dta$variable=="Took antibiotics",] %>%
    group_by(group) %>%
    do(test = fisher.test(matrix(c(.$n_shared[.$variable == "No antibiotics"], 
                                   .$n_shared[.$variable == "Took antibiotics"], 
                                   .$n_not_shared[.$variable == "No antibiotics"], 
                                   .$n_not_shared[.$variable == "Took antibiotics"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    )) 
  
  plot_dta <- filter(fig_dta, variable=="No antibiotics" | variable=="Took antibiotics")
  plot_dta <- merge(plot_dta, test_result, by=c("group"))
  plot_dta$variable <- factor(plot_dta$variable, levels=c("Took antibiotics", "No antibiotics"), labels=c("Yes", "No"))
  
  p_antib_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_antib_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_antib <- ggarrange(p_antib_1 + theme(plot.margin=unit(c(0.4, 0, 0, 0),"cm")), 
                       p_antib_2 + theme(plot.margin=unit(c(-0.2,0, 0.4, 0),"cm")), 
                       nrow=2, ncol=1, heights = c(1, 0.9))
  p_antib <- annotate_figure(p_antib, right=text_grob("Prenatal Antibiotic Use",  size = 8, rot = -90))
  
  
  # by delivery mode 
  test_result <- asv_dta[asv_dta$variable=="C-section" | asv_dta$variable=="Vaginal",] %>%
    group_by(group) %>%
    do(test = fisher.test(matrix(c(.$n_shared[.$variable == "Vaginal"], 
                                   .$n_shared[.$variable == "C-section"], 
                                   .$n_not_shared[.$variable == "Vaginal"], 
                                   .$n_not_shared[.$variable == "C-section"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))  
  plot_dta <- filter(fig_dta, variable=="C-section" | variable=="Vaginal")
  plot_dta <- merge(plot_dta, test_result, by=c("group"))
  
  p_delivery_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_delivery_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_delivery <- ggarrange(p_delivery_1 + theme(plot.margin=unit(c(0.4, 0, 0, 0),"cm")), 
                          p_delivery_2 + theme(plot.margin=unit(c(-0.2,0, 0.4, 0),"cm")), 
                          nrow=2, ncol=1, heights = c(1, 0.9))
  p_delivery <- annotate_figure(p_delivery, right=text_grob("Delivery Mode",  size = 8, rot = -90))
  
  
  
  # by ever breastfed
  test_result <- asv_dta[asv_dta$variable=="Yes" | asv_dta$variable=="No",] %>%
    group_by(group) %>%
    do(test = fisher.test(matrix(c(.$n_shared[.$variable == "No"], 
                                   .$n_shared[.$variable == "Yes"], 
                                   .$n_not_shared[.$variable == "No"], 
                                   .$n_not_shared[.$variable == "Yes"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  plot_dta <- filter(fig_dta, variable=="Yes" | variable=="No")
  plot_dta <- merge(plot_dta, test_result, by=c("group"))
  
  p_bf_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_bf_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_bf <- ggarrange(p_bf_1 + theme(plot.margin=unit(c(0.4, 0, 0, 0),"cm")), 
                    p_bf_2 + theme(plot.margin=unit(c(-0.2,0, 0.4, 0),"cm")), 
                    nrow=2, ncol=1, heights = c(1, 0.9))
  p_bf <- annotate_figure(p_bf, right=text_grob("Ever Breastfed", size = 8, rot = -90))
  
  
  p<- 
    ggarrange(ggarrange(p_overall, NULL, p_antib, widths = c(0.88, 0.15, 0.9), nrow=1), 
              ggarrange(p_delivery, NULL, p_bf, widths = c(1, 0.03, 0.9), nrow=1), nrow=2, heights = c(1, 1))
  
}


####  Version 2 for Aim 3
Bifido_Bac_genus_2 <- function(asv_dta){
  asv_dta <- asv_dta %>%
    group_by(cohort, name, variable) %>%
    mutate(sum_n=sum(n_pair))
  
  asv_dta <- asv_dta[asv_dta$type=="Shared",] %>% 
    group_by(group, variable) %>%
    summarise(n_detectable=sum(sum_n>0), n_shared=sum(n_pair>0), n_not_shared=sum(sum_n>0) - sum(n_pair>0))
  
  uplimt <- max(asv_dta$n_detectable) + 1
  
  fig_dta <- gather(asv_dta, type, n_pair, n_shared:n_not_shared, factor_key=TRUE)
  fig_dta$prop <- fig_dta$n_pair / fig_dta$n_detectable * 100
  fig_dta$type <- factor(fig_dta$type, levels=c("n_not_shared", "n_shared"), labels=c("Not shared", "Shared"))
  
  
  # overall
  plot_dta <- filter(fig_dta, variable=="Overall")
  
  p_overall_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota", ] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text = element_blank(),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) 
  
  p_overall_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota", ] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text = element_blank(),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank())
  
  
  p_overall <- ggarrange(p_overall_1 + theme(plot.margin=unit(c(0.5, 0, 0, 0),"cm")), 
                         p_overall_2 + theme(plot.margin=unit(c(-0.2,0, 0.5, 0),"cm")), 
                         nrow=2, ncol=1, heights = c(1, 0.9))
  p_overall <- annotate_figure(p_overall, right=text_grob("Overall",  size = 8, rot = -90))
  
  
  # by prenatal antibiotic use
  test_result <- asv_dta[asv_dta$variable=="No antibiotics" | asv_dta$variable=="Took antibiotics",] %>%
    group_by(group) %>%
    do(test = fisher.test(matrix(c(.$n_shared[.$variable == "No antibiotics"], 
                                   .$n_shared[.$variable == "Took antibiotics"], 
                                   .$n_not_shared[.$variable == "No antibiotics"], 
                                   .$n_not_shared[.$variable == "Took antibiotics"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    )) 
  
  plot_dta <- filter(fig_dta, variable=="No antibiotics" | variable=="Took antibiotics")
  plot_dta <- merge(plot_dta, test_result, by=c("group"))
  plot_dta$variable <- factor(plot_dta$variable, levels=c("Took antibiotics", "No antibiotics"), labels=c("Yes", "No"))
  
  p_antib_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_antib_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_antib <- ggarrange(p_antib_1 + theme(plot.margin=unit(c(0.4, 0, 0, 0),"cm")), 
                       p_antib_2 + theme(plot.margin=unit(c(-0.2,0, 0.4, 0),"cm")), 
                       nrow=2, ncol=1, heights = c(1, 0.9))
  p_antib <- annotate_figure(p_antib, right=text_grob("Prenatal Antibiotic Use",  size = 8, rot = -90))
  
  
  # by delivery mode 
  test_result <- asv_dta[asv_dta$variable=="C-section" | asv_dta$variable=="Vaginal",] %>%
    group_by(group) %>%
    do(test = fisher.test(matrix(c(.$n_shared[.$variable == "Vaginal"], 
                                   .$n_shared[.$variable == "C-section"], 
                                   .$n_not_shared[.$variable == "Vaginal"], 
                                   .$n_not_shared[.$variable == "C-section"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))  
  plot_dta <- filter(fig_dta, variable=="C-section" | variable=="Vaginal")
  plot_dta <- merge(plot_dta, test_result, by=c("group"))
  
  p_delivery_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_delivery_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_delivery <- ggarrange(p_delivery_1 + theme(plot.margin=unit(c(0.4, 0, 0, 0),"cm")), 
                          p_delivery_2 + theme(plot.margin=unit(c(-0.2,0, 0.4, 0),"cm")), 
                          nrow=2, ncol=1, heights = c(1, 0.9))
  p_delivery <- annotate_figure(p_delivery, right=text_grob("Delivery Mode",  size = 8, rot = -90))
  
  
  
  # by ever breastfed
  test_result <- asv_dta[asv_dta$variable=="Yes" | asv_dta$variable=="No",] %>%
    group_by(group) %>%
    do(test = fisher.test(matrix(c(.$n_shared[.$variable == "No"], 
                                   .$n_shared[.$variable == "Yes"], 
                                   .$n_not_shared[.$variable == "No"], 
                                   .$n_not_shared[.$variable == "Yes"]), 
                                 nrow = 2))) %>%
    ungroup() %>%
    mutate(p_value = map_dbl(test, "p.value"))
  
  # convert p-values to star symbols
  test_result <- test_result %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ "+",
      TRUE ~ ""
    ))
  
  plot_dta <- filter(fig_dta, variable=="Yes" | variable=="No")
  plot_dta <- merge(plot_dta, test_result, by=c("group"))
  
  p_bf_1 <- 
    plot_dta[plot_dta$group=="Sharing with maternal fecal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#0072B2")) + scale_y_discrete(position = "right") +
    theme_bw()+scale_x_reverse(limits=c(uplimt, 0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_bf_2 <- 
    plot_dta[plot_dta$group=="Sharing with maternal vaginal microbiota",] %>%
    ggplot(aes(y=variable, x=n_pair)) +
    geom_bar(aes(fill=type), stat="identity", position="stack", color="black", size=0.2) +  
    geom_text(aes(label=paste0(n_pair, " (", round(prop, 1), "%)"),  group=type), size = 2, position = position_stack(vjust=0.5))  +
    scale_fill_manual(values=c("white", "#E69F00")) + scale_y_discrete(position = "right") +
    theme_bw()+ scale_x_reverse(limits=c(uplimt, 0), position="top") + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),  axis.text.x = element_blank(), axis.text.y = element_text(size=6),
          legend.position="none", 
          panel.grid = element_blank(), panel.border = element_blank()) +
    geom_text(aes(x=0, y=1.5, label=significance), size=4, color="red", hjust=-0.2) 
  
  p_bf <- ggarrange(p_bf_1 + theme(plot.margin=unit(c(0.4, 0, 0, 0),"cm")), 
                    p_bf_2 + theme(plot.margin=unit(c(-0.2,0, 0.4, 0),"cm")), 
                    nrow=2, ncol=1, heights = c(1, 0.9))
  p_bf <- annotate_figure(p_bf, right=text_grob("Ever Breastfed", size = 8, rot = -90))
  
  
  p<- 
    ggarrange(ggarrange(p_overall, NULL, p_antib, widths = c(0.88, 0.15, 0.9), nrow=1), 
              ggarrange(p_delivery, NULL, p_bf, widths = c(1, 0.03, 0.9), nrow=1), nrow=2, heights = c(1, 1))
  
}



####=================   For Aim 3: Summarize number of shared ASV for each infant in the other four cohort  =====================####
aim3_asv_sharing_byinfant <- function(asv_table, i_id, m_id){
  ## Modify asv table
  asv_table <- as.data.frame(t(asv_table)) 
      # rows are samples, columns are taxa
  colnames(asv_table) <- asv_table[1, ]
      # assign colnames as ASV id
  asv_table <- asv_table[-(1:2), ] # remove unnecessary rows
  
  # only keep paired samples
  asv_table_paired <- asv_table %>% filter(row.names(asv_table) %in% as.list(i_id$sample.id) | row.names(asv_table) %in% as.list(m_id$sample.id))
  row_names <- row.names(asv_table_paired)
  
  asv_table_paired <- as.data.frame(lapply(asv_table_paired, as.integer)) 
    # transform counts to integer
  
  row.names(asv_table_paired) <- row_names
  
  
  ## Create ASV table separately for infant and mother
  asv_infant <- asv_table_paired %>% filter(row.names(asv_table_paired) %in% as.list(i_id$sample.id))
  asv_infant <- asv_infant[order(row.names(asv_infant)), ]
  
  asv_mother <- asv_table_paired %>% filter(row.names(asv_table_paired) %in% as.list(m_id$sample.id))
  asv_mother <- asv_mother[order(row.names(asv_mother)), ]
  
  # Create another relative abundance ASV table
  asv_infant_rlt <- asv_infant / rowSums(asv_infant) * 100
  asv_mother_rlt <- asv_mother / rowSums(asv_mother) * 100
  
  
  
  ## For each infant, find the matching mother
  matching_mother_indices <- sapply(rownames(asv_infant), function(infant_name) {
    infant_substr <- substr(infant_name, 7, 15)  # Extract characters 7 to 15 from the infant name
    matching_mother_index <- match(infant_substr, substr(rownames(asv_mother), 7, 15))
    return(matching_mother_index)
  })
  
  
  ## Create two matrices to store ASV sharing results
  share_res <- matrix(0, nrow = nrow(asv_infant), ncol = ncol(asv_infant)) # default 0 means not shared
  share_res_rlt_0.01 <- matrix(0, nrow = nrow(asv_infant), ncol = ncol(asv_infant)) # default 0 means not shared
  
  
  ## Loop through infants with matching mothers to determine the sharing for each ASV
  matching_rows <- !is.na(matching_mother_indices)
  for (n in 1:sum(matching_rows)) {
    infant_index <- n
    mother_index <- matching_mother_indices[infant_index]
    
    ## Definition 1: When present is defined as count>=1
    infant_vals <- asv_infant[infant_index, ]
    mother_vals <- asv_mother[mother_index, ]
    overlap_indices <- (infant_vals > 0) & (mother_vals > 0) #shared with mother
    notexist_indices <- (infant_vals == 0 & (mother_vals == 0)) # not existed
    
    ## Definition 4: When present is defined as relative abundance >=0.01%
    infant_vals_2 <- asv_infant_rlt[infant_index, ]
    mother_vals_2 <- asv_mother_rlt[mother_index, ]
    overlap_indices_2 <- (infant_vals_2 >= 0.01) & (mother_vals_2 >= 0.01) # shared with mother
    notexist_indices_2 <- (infant_vals_2 < 0.01) & (mother_vals_2 < 0.01) # not existed
    
    # Set the value to 3 when shared with mother
    share_res[infant_index, overlap_indices] <- 3 # shared 
    share_res[infant_index, notexist_indices] <- -1 # not existed
    
    share_res_rlt_0.01[infant_index, overlap_indices_2] <- 3 # shared
    share_res_rlt_0.01[infant_index, notexist_indices_2] <- -1 # not existed
  }
  
  ## Transfer to dataframe with row names being infant sample id
  share_res <- as.data.frame(share_res)
  rownames(share_res) <- rownames(asv_infant)
  
  share_res_rlt_0.01 <- as.data.frame(share_res_rlt_0.01)
  rownames(share_res_rlt_0.01) <- rownames(asv_infant_rlt)
  
  
  ##==========     Create datasets to summarize results by infant      ============##
  share_res_byinfant <- as.data.frame(i_id$sample.id)
  share_res_byinfant <- share_res_byinfant[order(share_res_byinfant$`i_id$sample.id`), ]
  share_res_byinfant <- as.data.frame(share_res_byinfant)
  colnames(share_res_byinfant) <- c("sample.id")
  
  
  share_res_byinfant <- share_res_byinfant %>%
    mutate(
      # count how many taxa are shared 
      shared = rowSums(share_res==3), 
      shared_0.01 = rowSums(share_res_rlt_0.01==3), 
      
      # count how many taxa existed in mother and infant
      sum_infant = rowSums(asv_infant>0),
      sum_infant_0.01 = rowSums(asv_infant_rlt>=0.01)
      
    ) 
  


  ##==========     Return results   =============##
  return(list(df1 = share_res_byinfant,  
              # this is summary table by infant
              df2 = share_res, df_3 = share_res_rlt_0.01)) 
              # these are results tables for each ASV and each infant, which will be used to create summary tables by ASV
}


##==========     For Aim3: Create datasets to summarize results by ASV      ============##
myfun <- function(asv_table, id.list, share_res, share_res_rlt_0.01){
  share_res <- filter(share_res, row.names(share_res) %in% as.list(id.list))
  share_res_rlt_0.01 <- filter(share_res_rlt_0.01, row.names(share_res_rlt_0.01) %in% as.list(id.list))
  
  share_res_byasv <- fun_asv_name(asv_table)
  
  share_res_byasv$shared = colSums(share_res==3)
  share_res_byasv$not_shared = colSums(share_res==0)
  
  share_res_byasv$shared_0.01 = colSums(share_res_rlt_0.01==3)
  share_res_byasv$not_shared_0.01 = colSums(share_res_rlt_0.01==0)
  
  return(share_res_byasv)
}

aim3_asv_sharing_byasv <- function(asv_table, i_id, result){
  asv_table <- as.data.frame(t(asv_table)) 
  # rows are samples, columns are taxa
  colnames(asv_table) <- asv_table[1, ]
  # assign colnames as ASV id
  asv_table <- asv_table[-(1:2), ] # remove unnecessary rows
  
  
  share_res <- result$df2 
  share_res_rlt_0.01 <- result$df_3
  
  
  # Overall
  id.list <- as.list(i_id$sample.id)
  share_res_byasv1 <- myfun(asv_table, id.list, share_res, share_res_rlt_0.01)
  share_res_byasv1$variable <- "Overall"

  # By prenatal antibiotics
  id.list <- as.list(i_id[i_id$med_antib_pn=="No antibiotics",]$sample.id)
  share_res_byasv4 <- myfun(asv_table, id.list, share_res, share_res_rlt_0.01)
  share_res_byasv4$variable <- "No"

  id.list <- as.list(i_id[i_id$med_antib_pn=="Took antibiotics",]$sample.id)
  share_res_byasv5 <- myfun(asv_table, id.list, share_res, share_res_rlt_0.01)
  share_res_byasv5$variable <- "Yes"

  # By delivery mode
  id.list <- as.list(i_id[i_id$birth_delivm=="Vaginal delivery",]$sample.id)
  share_res_byasv2 <- myfun(asv_table, id.list, share_res, share_res_rlt_0.01)
  share_res_byasv2$variable <- "Vaginal"

  id.list <- as.list(i_id[i_id$birth_delivm=="Cesarean delivery",]$sample.id)
  share_res_byasv3 <- myfun(asv_table, id.list, share_res, share_res_rlt_0.01)
  share_res_byasv3$variable <- "C-section"

  # Combine
  share_res_byasv <- rbind(share_res_byasv1, share_res_byasv4, share_res_byasv5, share_res_byasv2, share_res_byasv3)
  
  return(share_res_byasv)
}
