#code for creating graphs
##please change the paths accordingly

##step 1.
##distribution of many unique credible sets per locus (careful! it could be multiple ones) per SNP size

dat_new <- read_csv('/Users/koromm03/Downloads/snakemake_finemap_pgc3_bip/processing/sizeCS/polyfun_finemap_HRC_windows_sizeCS.csv')
dat2 <- dat_new[!(is.na(dat_new$locus)),]

sum(dat2$'1' ==1, na.rm = TRUE)
sum(dat2$'1' >= 2 & dat2$'1' <= 5, na.rm= TRUE)
sum(dat2$'1' >= 6 & dat2$'1' <= 10, na.rm= TRUE)
sum(dat2$'1' >= 11 & dat2$'1' <= 20, na.rm= TRUE)
sum(dat2$'1' >= 21 & dat2$'1' <= 50, na.rm= TRUE)
sum(dat2$'1' >= 50, na.rm=TRUE)
  
#repeat for CSs 2,3,4, etc
sum(dat2$'2' ==1, na.rm = TRUE)
sum(dat2$'2' >= 2 & dat2$'2' <= 5, na.rm= TRUE)
sum(dat2$'2' >= 6 & dat2$'2' <= 10, na.rm= TRUE)
sum(dat2$'2' >= 11 & dat2$'2' <= 20, na.rm= TRUE)
sum(dat2$'2' >= 21 & dat2$'2' <= 50, na.rm= TRUE)
sum(dat2$'2' >= 50, na.rm=TRUE)

sum(dat2$'3' ==1, na.rm = TRUE)
sum(dat2$'3' >= 2 & dat2$'3' <= 5, na.rm= TRUE)
sum(dat2$'3' >= 6 & dat2$'3' <= 10, na.rm= TRUE)
sum(dat2$'3' >= 11 & dat2$'3' <= 20, na.rm= TRUE)
sum(dat2$'3' >= 21 & dat2$'3' <= 50, na.rm= TRUE)
sum(dat2$'3' >= 50, na.rm=TRUE)

sum(dat2$'4' ==1, na.rm = TRUE)
sum(dat2$'4' >= 2 & dat2$'4' <= 5, na.rm= TRUE)
sum(dat2$'4' >= 6 & dat2$'4' <= 10, na.rm= TRUE)
sum(dat2$'4' >= 11 & dat2$'4' <= 20, na.rm= TRUE)
sum(dat2$'4' >= 21 & dat2$'4' <= 50, na.rm= TRUE)
sum(dat2$'4' >= 50, na.rm=TRUE)

sum1 = sum(dat2$'1' ==1, na.rm = TRUE) + sum(dat2$'2' ==1, na.rm = TRUE) + sum(dat2$'3' ==1, na.rm = TRUE) + sum(dat2$'4' ==1, na.rm = TRUE)
sum2 = sum(dat2$'1' >= 2 & dat2$'1' <= 5, na.rm= TRUE) + sum(dat2$'2' >= 2 & dat2$'2' <= 5, na.rm= TRUE) + sum(dat2$'3' >= 2 & dat2$'3' <= 5, na.rm= TRUE) + sum(dat2$'4' >= 2 & dat2$'4' <= 5, na.rm= TRUE)
sum3 = sum(dat2$'1' >= 6 & dat2$'1' <= 10, na.rm= TRUE)  + sum(dat2$'2' >= 6 & dat2$'2' <= 10, na.rm= TRUE) + sum(dat2$'3' >= 6 & dat2$'3' <= 10, na.rm= TRUE) + sum(dat2$'4' >= 6 & dat2$'4' <= 10, na.rm= TRUE) 
sum4 = sum(dat2$'1' >= 11 & dat2$'1' <= 20, na.rm= TRUE) + sum(dat2$'2' >= 11 & dat2$'2' <= 20, na.rm= TRUE) + sum(dat2$'3' >= 11 & dat2$'3' <= 20, na.rm= TRUE) + sum(dat2$'4' >= 11 & dat2$'4' <= 20, na.rm= TRUE)
sum5 = sum(dat2$'1' >= 21 & dat2$'1' <= 50, na.rm= TRUE) + sum(dat2$'2' >= 21 & dat2$'2' <= 50, na.rm= TRUE) + sum(dat2$'3' >= 21 & dat2$'3' <= 50, na.rm= TRUE) + sum(dat2$'4' >= 21 & dat2$'4' <= 50, na.rm= TRUE)
sum6 = sum(dat2$'1' >= 50, na.rm=TRUE) + sum(dat2$'2' >= 50, na.rm=TRUE) + sum(dat2$'3' >= 50, na.rm=TRUE) + sum(dat2$'4' >= 50, na.rm=TRUE)

#create a dataset
df <- tibble(SNPs_in_CS = c("size_1", "size_2-5", "size_6-10", "size_11-20", "size_21-50", "size_50+"),
               finemap_loci = c(sum1, sum2, sum3, sum4, sum5, sum6))  #poly-susie-ukb-windows
#df <- tibble(SNPs_in_CS = c("size_1", "size_2-5", "size_6-10", "size_11-20", "size_21-50", "size_50+"),
             #finemap_loci = c(5,2,7,11,15,11))  #poly_susie-hrc-windows  

#to reorder stacks, check library forcats and do: fill = fct_reorder
library(forcats)
df$SNPs_in_CS <- fct_reorder(df$SNPs_in_CS, df$finemap_loci)
group.colors <- c("size_1"= "#F8766D", "size_2-5" = "#B79F00", "size_6-10"= "#00BA38" , "size_11-20"= "#00BFC4", "size_21-50"= "#619CFF", "size_50+"= "#F564E3")

#Grouped
 c <- ggplot(df, aes(x = "", y = finemap_loci, fill = SNPs_in_CS)) +
    geom_col(width=0.2) +
    geom_text(aes(label = paste0(finemap_loci)), 
              position = position_stack(vjust = 0.5)) + labs(x="", y= "N of finemapped loci") +
   scale_fill_manual(values=group.colors) + coord_flip()
   
   
 ##step 2.
 ##code borrowed by Ashvin Ravi (Raj Lab)
 ##smallest CS per locus for all fine-mapping methods given a certain LD panel and windows range

library(ggplot2)
library(dplyr)
library(tidyverse)

concatenate_results <- function(fldr, method) {
  files = list.files(fldr, pattern = '.gz')
  new_results = data.frame()
  for (f in files) {
    setwd(fldr)
    file = read_tsv(f)
    file$LOCUS = gsub(".*\\.(.*)\\..*", "\\1", f)
    print(gsub(".*\\.(.*)\\..*", "\\1", f))
    new_results = rbind(new_results, file)
  }
  return(new_results)
}

finemap_results = concatenate_results('/sc/arion/projects/ad-omics/ashvin/finemapping/polyfun/ASHG_finemapping_results/finemap', finemap)
polyfun_susie_results = concatenate_results('/sc/arion/projects/ad-omics/ashvin/finemapping/polyfun/ASHG_finemapping_results/polyfun_susie/', polyfun_susie)
polyfun_finemap_results = concatenate_results('/sc/arion/projects/ad-omics/ashvin/finemapping/polyfun/ASHG_finemapping_results/polyfun_finemap', polyfun_finemap)
susie_results = concatenate_results('/sc/arion/projects/ad-omics/ashvin/finemapping/polyfun/ASHG_finemapping_results/susie/')
susie_ld_gwas = concatenate_results('/sc/arion/projects/ad-omics/ashvin/finemapping/polyfun/ASHG_finemapping_results/susie_gwas_ld/')

credible_set_graph <- function(dataset, method, threshold) {
  credible_set_finemap_results <- dataset %>% group_by(LOCUS, CREDIBLE_SET) %>% dplyr::filter(PIP >= threshold) %>% summarise(sum(CREDIBLE_SET))
  credible_set_finemap_results$xaxis <- method
  print(head(credible_set_finemap_results))
  if (length(dataset$LOCUS[!dataset$LOCUS %in% credible_set_finemap_results$LOCUS]) > 0) {
    credible_set_zero_under_threshold <- data.frame(dataset$LOCUS[!dataset$LOCUS %in% credible_set_finemap_results$LOCUS])
    dim(credible_set_zero_under_threshold)
    credible_set_zero_under_threshold$credible_set_sum = 0
    print(head(credible_set_zero_under_threshold))
    colnames(credible_set_zero_under_threshold) = c('LOCUS', 'credible_set_sum')
    credible_set_zero_under_threshold <- unique(credible_set_zero_under_threshold)
    credible_set_finemap_results_final <- credible_set_finemap_results %>% dplyr::filter(CREDIBLE_SET != 0) %>% summarize(min(`sum(CREDIBLE_SET)`)) 
    zero_cs <- credible_set_finemap_results[!credible_set_finemap_results$LOCUS %in% credible_set_finemap_results_final$LOCUS,]
    zero_cs <- zero_cs[, c(1,3)]
    colnames(credible_set_finemap_results_final) = c('LOCUS', 'credible_set_sum')
    colnames(zero_cs) = c('LOCUS', 'credible_set_sum')
    credible_set_finemap_results_final <- rbind(zero_cs, credible_set_zero_under_threshold, credible_set_finemap_results_final)
  } else {
    credible_set_finemap_results_final <- credible_set_finemap_results %>% dplyr::filter(CREDIBLE_SET != 0) %>% summarize(min(`sum(CREDIBLE_SET)`)) 
    zero_cs <- credible_set_finemap_results[!credible_set_finemap_results$LOCUS %in% credible_set_finemap_results_final$LOCUS,]
    zero_cs <- zero_cs[, c(1,3)]
    colnames(credible_set_finemap_results_final) = c('LOCUS', 'credible_set_sum')
    colnames(zero_cs) = c('LOCUS', 'credible_set_sum')
    credible_set_finemap_results_final <- rbind(zero_cs, credible_set_finemap_results_final)
  }
  # =credible_set_finemap_results <- rbind(credible_set_finemap_results_final, zero_cs)
  
  credible_set_finemap_results_final$category = '15+'
  credible_set_finemap_results_final$decoy2 <- 9
  
  
  credible_set_finemap_results_final$category[credible_set_finemap_results_final$credible_set_sum < 15 & credible_set_finemap_results_final$credible_set_sum >= 11] = '11-14'
  credible_set_finemap_results_final$decoy2[credible_set_finemap_results_final$credible_set_sum < 15 & credible_set_finemap_results_final$credible_set_sum >= 11] = 8
  
  credible_set_finemap_results_final$category[credible_set_finemap_results_final$credible_set_sum <= 10 & credible_set_finemap_results_final$credible_set_sum >= 8] = '8-10'
  credible_set_finemap_results_final$decoy2[credible_set_finemap_results_final$credible_set_sum <= 10 & credible_set_finemap_results_final$credible_set_sum >= 8] = 6
  
  credible_set_finemap_results_final$category[credible_set_finemap_results_final$credible_set_sum <= 7 & credible_set_finemap_results_final$credible_set_sum >= 5] = '5-7'
  credible_set_finemap_results_final$decoy2[credible_set_finemap_results_final$credible_set_sum <= 7 & credible_set_finemap_results_final$credible_set_sum >= 5] = 5
  
  credible_set_finemap_results_final$category[credible_set_finemap_results_final$credible_set_sum <= 4 & credible_set_finemap_results_final$credible_set_sum >= 2] = '2-4'
  credible_set_finemap_results_final$decoy2[credible_set_finemap_results_final$credible_set_sum <= 4 & credible_set_finemap_results_final$credible_set_sum >= 2] = 3
  
  credible_set_finemap_results_final$category[credible_set_finemap_results_final$credible_set_sum == 1] = '1'
  credible_set_finemap_results_final$decoy2[credible_set_finemap_results_final$credible_set_sum == 1] = 2
  
  credible_set_finemap_results_final$category[credible_set_finemap_results_final$credible_set_sum == 0] = '0'
  credible_set_finemap_results_final$decoy2[credible_set_finemap_results_final$credible_set_sum == 0] = 0
  print(head(credible_set_finemap_results_final[credible_set_finemap_results_final$category == '15+',]))
  credible_set_finemap_results_final$xaxis <- 1
  credible_set_finemap_results_final$decoy <- seq(1, nrow(credible_set_finemap_results_final), by=1)
  
  
  
  temp <- data.frame(table(credible_set_finemap_results_final$category))
  print(head(temp))
  
  temp$Var1 <- factor(temp$Var1, 
                      levels = c('15+', 
                                 '11-14',
                                 '8-10',
                                 '5-7',
                                 '2-4',
                                 '1',
                                 '0'))
  
  sp <- ggplot(data=temp, 
               aes(x=method, 
                   y=Freq, 
                   fill=Var1,
                   label=Freq
               )
  ) + 
    geom_bar(position='stack', 
             stat='identity') + 
    geom_text(aes(method, label = Freq), data = temp, position=position_stack(vjust=0.5)) +
    coord_flip() + 
    ggtitle(method) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.05), legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(5,5,5,5)) +
    scale_fill_manual(values=c("#D53E4F", 
                                        "#FC8D59", 
                                        "#FEE08B", 
                                        "#FFFFBF",
                                        "#E6F598", 
                                        "#99D594", 
                                        "#3288BD"),
                                        drop=FALSE, 
                      name='CS size', 
                      labels=c("15+", 
                               "11-14", 
                               "8-10", 
                               "5-7", 
                               "2-4", 
                               "1", 
                               "0")) 
  return(sp)
}

finemap_credible_sets <- credible_set_graph(finemap_results, 'FINEMAP', 0)
polyfun_finemap_credible_sets <- credible_set_graph(polyfun_finemap_results, 'PolyFun + FINEMAP', 0)
susie_credible_sets <- credible_set_graph(susie_results, 'SuSiE', 0)
polyfun_susie_credible_sets <- credible_set_graph(polyfun_susie_results, 'PolyFun + SuSiE', 0)
susie_ld_gwas_credible_sets <- credible_set_graph(susie_ld_gwas, 'SuSiE + GWAS LD', 0)

CS_graph_zero <- ggpubr::ggarrange(finemap_credible_sets, polyfun_finemap_credible_sets, susie_credible_sets, polyfun_susie_credible_sets, susie_ld_gwas_credible_sets, ncol=1, nrow=5, common.legend = TRUE, legend="top")

finemap_credible_sets_0.001 <- credible_set_graph(finemap_results, 'FINEMAP', 0.001)
polyfun_finemap_credible_sets_0.001 <- credible_set_graph(polyfun_finemap_results, 'PolyFun + FINEMAP', 0.001)
susie_credible_sets_0.001 <- credible_set_graph(susie_results, 'SuSiE', 0.001)
polyfun_susie_credible_sets_0.001 <- credible_set_graph(polyfun_susie_results, 'PolyFun + SuSiE', 0.001)
susie_ld_gwas_credible_sets_0.001 <- credible_set_graph(susie_ld_gwas, 'SuSiE + GWAS LD', 0.001)

CS_graph_0.001 <- ggpubr::ggarrange(finemap_credible_sets_0.001, polyfun_finemap_credible_sets_0.001, susie_credible_sets_0.001, polyfun_susie_credible_sets_0.001, susie_ld_gwas_credible_sets_0.001, ncol=1, nrow=5, common.legend = TRUE, legend="right")

finemap_credible_sets_0.1 <- credible_set_graph(finemap_results, 'FINEMAP', 0.1)
polyfun_finemap_credible_sets_0.1 <- credible_set_graph(polyfun_finemap_results, 'PolyFun + FINEMAP', 0.1)
susie_credible_sets_0.1 <- credible_set_graph(susie_results, 'SuSiE', 0.1)
polyfun_susie_credible_sets_0.1 <- credible_set_graph(polyfun_susie_results, 'PolyFun + SuSiE', 0.1)
susie_ld_gwas_credible_sets_0.1 <- credible_set_graph(susie_ld_gwas, 'SuSiE + GWAS LD', 0.1)

CS_graph_0.1 <- ggpubr::ggarrange(finemap_credible_sets_0.1, polyfun_finemap_credible_sets_0.1, susie_credible_sets_0.1, polyfun_susie_credible_sets_0.1, susie_ld_gwas_credible_sets_0.1, ncol=1, nrow=5, common.legend = TRUE, legend="right")


##step 3.
##heatmap as in Figure 3 of the fine-mapping paper

dat <- read_csv('~/Desktop/scripts/visualisation/heatmap_apr23.csv')
dat$variables = dat$`Genetic ID`
dat$`Genetic ID` = NULL
head(dat)
test <- dat %>% pivot_longer(-variables,names_to = "rsIDs", values_to = "Value")

test$variables[test$variables=="Splicing QTLs (ROSMAP)"] <- "Brain sQTLs (ROSMAP)"
test$variables[test$variables=="Splicing QTLs (CMC)"] <-"Brain sQTLs (CMC)"
test$variables[test$variables=="BipEX"] <-"BipEx"

# ggplot(test, aes(variables, rsIDs)) +
#   geom_tile(aes(fill = Value), colour = "grey50") + theme(axis.text.y = element_text(size = 9)) +
#   theme(axis.text.x = element_text(size = 9, angle = 35, hjust = 1)) + ylab("Union Consensus SNPs") +xlab("") +
#   theme(legend.title=element_blank()) + scale_fill_gradient2(na.value= "yellow")

p <- test %>% mutate(variables = fct_relevel(variables, "Astrocyte Enhancers", "Astrocyte Promoters", "Neuronal Enhancers", "Neuronal Promoters", "Brain eQTLs (PEC)", "Brain sQTLs (CMC)", "Brain sQTLs (ROSMAP)", "BipEx", "pLI")) %>%
  ggplot(aes(variables, rsIDs)) +
  geom_tile(aes(fill = Value), colour = "grey50") + theme(axis.text.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 9, angle = 25, hjust = 1)) + ylab("Union Consensus SNPs") + xlab("") +
  theme(legend.title=element_blank()) + scale_fill_gradient2(na.value= "grey70")
  
  
##step 4.
##locus plots ~ case example for FURIN

library(readr)
library(ggbio)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(LDlinkR)
library(valr)

        
plot_variants_in_ld <- function(ldlink_results, start, end) {
  ldlink_sp <- ggplot(ldlink_results, aes(x=BP, y=-log10(PVALUE))) +
    geom_point(aes(colour = R2), binaxis='y', stackdir='center') +
    scale_colour_gradientn(colours = c("blue","white","red")) +
    xlab(element_blank()) + ylab('-log10(GWAS p-value)') +
    scale_x_continuous(breaks = seq(start, end, by = 50000)) +
    coord_cartesian(xlim = c(start, end)) +
    theme_bw()
}

plot_fine_mapping_variants <- function(finemapping_results, start, end) {
  finemapping_results$CREDIBLE_SET <- as.character(finemapping_results$CREDIBLE_SET)
  a <- ggplot(finemapping_results, aes(x=start, y=PIP), ylim = c(0,1)) +
    ggplot2::geom_bar(aes(color = CREDIBLE_SET, fill = CREDIBLE_SET), stat = 'identity', binaxis='y', stackdir='center') +
    xlab(element_blank()) +
    scale_fill_discrete(drop=TRUE, limits = levels(finemapping_results$CREDIBLE_SET)) +
    scale_x_continuous(breaks = seq(start, end, by = 50000)) +
    scale_y_continuous(breaks = seq(0, 100, by = 0.25), limits = c(0, 1)) +
    coord_cartesian(xlim = c(start, end)) +
    theme_bw()

  return(a)
}

nott_epigenome_plot <- function(snp_loc, start, end, color){
  plot <- ggplot() + geom_rect(snp_loc, mapping=aes(xmin=start,xmax=end, ymin=0, ymax=1, fill = label)) +
    scale_fill_manual(values = c(color)) +
    #guides(fill = FALSE) +
    scale_x_continuous(breaks = seq(start, end, by = 50000)) +
    coord_cartesian(xlim = c(start, end)) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  return(plot)
}

abc_model_plot <- function(plac_loc, ABC_start, ABC_end, TSS_start) {

  ## make plot
  ggplot() + ggbio::geom_arch(data = plac_loc, aes(x = end1, xend = start2, colour = overlap, alpha = overlap) ) + theme_classic() + labs(x = "", y = "") +
    ggbio::geom_rect(data = plac_loc, aes(xmin = start1, xmax = end1, fill = overlap, alpha = overlap ), ymin = 0, ymax = 1, colour = NA ) +
    ggbio::geom_rect(data = plac_loc, aes(xmin = start2, xmax = start2 + 1000, fill = overlap, alpha = overlap), ymin = 0, ymax = 1, colour = NA) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank() ) +
    scale_fill_manual(values = c("darkgray", "black")) +
    scale_colour_manual(values = c("gray", "black")) +
    scale_alpha_manual(values = c(0.1, 1)) +
    guides(colour = FALSE, alpha = FALSE, fill = FALSE)
}


##ensembldb

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 75)

db.gr <- ensembldb::transcripts(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) %>%
  data.table::as.data.table() %>%
  dplyr::filter( tx_biotype == "protein_coding") %>%
  #dplyr::mutate(index=row.names(.)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::slice_max(width, n = 1)


gene_plot <- function(chr_decoy, start_decoy, end_decoy){
  
  snp_loc_start <- start_decoy
  snp_loc_end <- end_decoy
  snp_gr <- GenomicRanges::GRanges(seqnames = chr_decoy, ranges = IRanges::IRanges(start = snp_loc_start, end = snp_loc_end))
  GenomeInfoDb::seqlevelsStyle(snp_gr) <- "NCBI"
  
  db_loc <- subset(db.gr,
                   seqnames == chr_decoy &
                     (start >= snp_loc_start | end <= snp_loc_end + 4000000) ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(db_loc) <- "NCBI"
  
  # db_loc$symbol <- factor(db_loc$symbol, levels = unique(db_loc$symbol), ordered = T)
  edb <-  ensembldb::addFilter( EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                AnnotationFilter::TxIdFilter(db_loc$tx_id) )
  
  plot <- ggplot() + ggbio::geom_alignment(edb, which = snp_gr, rect.height = unit(3, "cm"),
                                           names.expr = "gene_name",
                                           aes(group = gene_name,
                                               fill='arrow',
                                               arrow.rate = 1,
                                               length = unit(0.001, "cm"),
                                               show.legend=FALSE)) + scale_fill_manual(values = c("#7A7A7A")) +
    theme_classic() + theme(axis.text = element_text(colour = "black", size = 14)) + scale_x_continuous(breaks = seq(start_decoy, end_decoy, by = 50000)) +
    coord_cartesian(xlim = c(start_decoy, end_decoy))
  
  return(plot)
}

# FURIN GENE

bip_gwas <- read_tsv('/Users/koromm03/Desktop/locus_plots/Daner_2020.processed.tsv.gz')
FURIN_variants_in_ld <- LDproxy('rs4702', pop = "EUR", r2d = "r2", token = '1f09218e3faa')

FURIN_variants_in_ld <- FURIN_variants_in_ld[, c(1,7)]
names(FURIN_variants_in_ld)[names(FURIN_variants_in_ld) == 'RS_Number'] = 'SNP'
FURIN_variants_to_plot <- inner_join(FURIN_variants_in_ld, bip_gwas)
names(FURIN_variants_to_plot)[names(FURIN_variants_to_plot) == 'P'] = 'PVALUE'
FURIN_gwas_plot <- plot_variants_in_ld(FURIN_variants_to_plot, 91375000, 91475000)
#add +geom_label or geom_point to hihglight SNPs above a certain Pval
# FURIN_gwas_plot + geom_point()  + geom_text(data=subset(FURIN_variants_to_plot, PVALUE < 1e-06))
#add GWAS line

bip_finemapping_results <- read_csv('/Users/koromm03/Desktop/polyfun_processing/merged_dfs/only_susie_HRC.csv')
names(bip_finemapping_results)[names(bip_finemapping_results) == 'BP'] = 'start'
names(bip_finemapping_results)[names(bip_finemapping_results) == 'CHR'] = 'chrom'
bip_finemapping_results$end <- bip_finemapping_results$start + 1
bip_finemapping_results$chrom <- paste('chr', bip_finemapping_results$chrom, sep = '')
FURIN_finemapping_results <- bip_finemapping_results %>% dplyr::filter(LOCUS == '91426560')
names(FURIN_finemapping_results)[names(FURIN_finemapping_results) == 'BP'] = 'start'
FURIN_finemap_plot <- plot_fine_mapping_variants(FURIN_finemapping_results, 91375000, 91475000)
#add +geom_label or geom_point to highlight SNPs above a certain PIP

bip_finemapping_results2 <- read_csv('/Users/koromm03/Desktop/polyfun_processing/merged_dfs/polyfun_susie_HRC.csv')
names(bip_finemapping_results2)[names(bip_finemapping_results2) == 'BP'] = 'start'
names(bip_finemapping_results2)[names(bip_finemapping_results2) == 'CHR'] = 'chrom'
bip_finemapping_results2$end <- bip_finemapping_results2$start + 1
bip_finemapping_results2$chrom <- paste('chr', bip_finemapping_results2$chrom, sep = '')
FURIN_finemapping_results2 <- bip_finemapping_results2 %>% dplyr::filter(LOCUS == '91426560')
names(FURIN_finemapping_results2)[names(FURIN_finemapping_results2) == 'BP'] = 'start'
FURIN_finemap_plot2 <- plot_fine_mapping_variants(FURIN_finemapping_results2, 91375000, 91475000)

bip_finemapping_results3 <- read_csv('/Users/koromm03/Desktop/polyfun_processing/merged_dfs/only_finemap_HRC.csv')
names(bip_finemapping_results3)[names(bip_finemapping_results3) == 'BP'] = 'start'
names(bip_finemapping_results3)[names(bip_finemapping_results3) == 'CHR'] = 'chrom'
bip_finemapping_results3$end <- bip_finemapping_results3$start + 1
bip_finemapping_results3$chrom <- paste('chr', bip_finemapping_results3$chrom, sep = '')
FURIN_finemapping_results3 <- bip_finemapping_results3 %>% dplyr::filter(LOCUS == '91426560')
names(FURIN_finemapping_results3)[names(FURIN_finemapping_results3) == 'BP'] = 'start'
FURIN_finemap_plot3 <- plot_fine_mapping_variants(FURIN_finemapping_results3, 91375000, 91475000)

bip_finemapping_results4 <- read_csv('/Users/koromm03/Desktop/polyfun_processing/merged_dfs/polyfun_finemap_HRC.csv')
names(bip_finemapping_results4)[names(bip_finemapping_results4) == 'BP'] = 'start'
names(bip_finemapping_results4)[names(bip_finemapping_results4) == 'CHR'] = 'chrom'
bip_finemapping_results4$end <- bip_finemapping_results4$start + 1
bip_finemapping_results4$chrom <- paste('chr', bip_finemapping_results4$chrom, sep = '')
FURIN_finemapping_results4 <- bip_finemapping_results4 %>% dplyr::filter(LOCUS == '91426560')
names(FURIN_finemapping_results4)[names(FURIN_finemapping_results4) == 'BP'] = 'start'
FURIN_finemap_plot4 <- plot_fine_mapping_variants(FURIN_finemapping_results4, 91375000, 91475000)


neuronal_promoter <- read.table('/Users/koromm03/Desktop/locus_plots/Neuronal_promoters.bed', sep = '\t', header = FALSE)
colnames(neuronal_promoter) <- c('chrom', 'start', 'end', 'peak')
FURIN_promoter_peaks <- neuronal_promoter %>% dplyr::filter(chrom == 'chr15', start >= 91375000)
FURIN_promoter_peaks <- FURIN_promoter_peaks %>% dplyr::filter(chrom == 'chr15', end <= 91475000)
FURIN_promoter_peaks$label <- 'promoter'
a = nott_epigenome_plot(FURIN_promoter_peaks, 91375000, 91475000, "#7F00FF")

neuronal_enhancer <- read.table('/Users/koromm03/Desktop/locus_plots/Neuronal_enhancers.bed', sep = '\t', header = FALSE)
colnames(neuronal_enhancer) <- c('chrom', 'start', 'end', 'peak')
FURIN_enhancer_peaks <- neuronal_enhancer %>% dplyr::filter(chrom == 'chr15', start >= 91375000)
FURIN_enhancer_peaks <- FURIN_enhancer_peaks %>% dplyr::filter(chrom == 'chr15', end <= 91475000)
FURIN_enhancer_peaks$label <- 'enhancer'
b = nott_epigenome_plot(FURIN_enhancer_peaks, 91375000, 91475000, "#DAA520")

abc_model_enhancers <- read_tsv('/Users/koromm03/Desktop/locus_plots/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz')
bipolar_iPSC <- abc_model_enhancers[abc_model_enhancers$CellType == 'bipolar_neuron_from_iPSC-ENCODE',]

names(bipolar_iPSC)[names(bipolar_iPSC) == 'chr'] = 'chrom'
bip_neuron_intersect <- bed_intersect(bipolar_iPSC, bip_finemapping_results)
bip_neuron_intersect <- bip_neuron_intersect[, c(1,2,3,7,8,21,24,26,27,34,39)]

FURIN_ABC_enhancers <- bipolar_iPSC %>% dplyr::filter(chrom == 'chr15',
                                                                      start >= 91375000,
                                                                      end <= 91475000)

FURIN_ABC_enhancers$falls_in <- 'no'
FURIN_ABC_enhancers$falls_in[FURIN_ABC_enhancers$start <= 91426560 & FURIN_ABC_enhancers$end >= 91426560] = 'yes'

FURIN_abc_graph <- ggplot() + ggbio::geom_arch(data = FURIN_ABC_enhancers,
                                               aes(x = (start + end)/2,
                                                   xend = TargetGeneTSS,
                                                   color = falls_in,
                                                   alpha = ABC.Score)) +
  theme_classic() + labs(x = "", y = "") +
  scale_fill_manual(values = c("darkgray", "black")) +
  scale_colour_manual(values = c("darkgray", "blue")) +
  ggbio::geom_rect(data = FURIN_ABC_enhancers, aes(xmin = start, xmax = end, alpha = ABC.Score), ymin = 0, ymax = 1, colour = NA ) +
  ggbio::geom_rect(data = FURIN_ABC_enhancers, aes(xmin = TargetGeneTSS, xmax = TargetGeneTSS + 1000, alpha = ABC.Score), ymin = 0, ymax = 1, fill = 'purple') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank() ) +
  scale_x_continuous(breaks = seq(91375000, 91475000, by = 50000)) +
  coord_cartesian(xlim = c(91375000, 91475000))

FURIN_plot <- gene_plot(15,91375000, 91475000)

FURIN <- cowplot::plot_grid(FURIN_gwas_plot, FURIN_finemap_plot, FURIN_finemap_plot2, FURIN_finemap_plot3, FURIN_finemap_plot4, a, b, FURIN_abc_graph, FURIN_plot, align='v', nrow=9, rel_heights = c(1,.5,.5,.5,.5,.5,.5,.6,1.2))




