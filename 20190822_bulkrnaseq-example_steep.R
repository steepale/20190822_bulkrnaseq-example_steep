#'---
#' title: "RNASeq Example Code: Exploratory data analysis reveals surrogate variables -- Don't get burned by the gonads"
#' author: Alec Steep
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         
#'         highlight: zenburn
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' # Example code from a bulk RNA Seq Analysis of Marek's disease lymphomas seeded in gonadal tissue
#' 
#' 
#' ## Explanation of Samples and Data 
#' ### Explaination of Biological Samples and Sequencing
#' Biological samples obtained from Marek's disease (MD) lymphomas seeded in the gonads of males and females of F1 birds--cross of parental lines 6 (resistent to MD) and 7 (suspectible to MD)--infected at hatch with 1000 pfu of Marek's disease virus (MDV) (strain JM102W). Birds were euthanized and tumors collected when birds were moribund or if they reached eight weeks of age (birds challenged at hatch). RNA was extracted via the miTNeasy mini kit (Qiagen). Sequencing was performed in 1 batch; 26 gonadal tumor samples from 22 tumors from 22 birds and 8 samples of isolated CD4+ spleenic T-cells from 8 uninfected birds on 05/28/2015. Experimental design decision:normal CD4+ T cell populations cannot be obtained from MD-infected birds; therefore, normal CD4+ cell populations were obtained from uninfected birds of 6x7 F1 cross. RNA-Seq libraries were prepared using the NuGen Ovation Single Cell RNA-Seq System (stranded). All samples underwent sequencing to produce 125bp paired-end reads via Illumina HiSeq machines at the Michigan State University RTSF Genomics Core. Project designed and performed by Hans Cheng (Hans.Cheng@ars.usda.gov) and Alec Steep (alec.steep@gmail.com).
#'
#'
#' ### Sample/File Nomenclature
#'
#' #### Primary Tumor Samples:
#'`017834-2` 
#'`BIRD#-TUMOR#`
#'`017834`: Bird ID  
#'`2`: Tumor ID relative to bird (e.g. tumor 2 collected from bird 017834)  
#'
#' #### Replicate Tumor Samples:
#'`017901-2_2`  
#'`BIRD#-TUMOR#_Replicate`
#'`017901`: Bird ID  
#'`2`: Tumor ID relative to bird (e.g. tumor 2 collected from bird 017901)  
#'`_2`: Identifier of biological replicate from distant side of tumor  
#'
#' #### Primary Control Samples:
#' `017834-0`
#' `BIRD#-NORMAL` (always 0)
#' `017834`: Bird ID  
#' `0`: Normal ID (always 0)
#'
#'
#' ## Setup the Environment

#+ Setup Environment

############################################################
##### Dependencies #########################################
############################################################

# Load dependencies
#################
#BiocManager::install("SeqGSEA", version = "3.8")
pacs...man <- c("GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2",
                "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr",
                "org.Gg.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt","feather",
                "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer')

# Install packages if need be (NTD: build docker image with dependencies)
new.packages <- pacs...man[!(pacs...man %in% installed.packages()[,"Package"])]
if( length(new.packages) ) 
        install.packages(new.packages)

lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})
#################

############################################################
##### Functions ############################################
############################################################

# Make the 'not in' operator
############################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
############################

# Function to relabel RNASeq read names to orthologs
###############################################################################################
chick2ortho_df <- function(x) {
        # Ensure 'x' is a data.frame
        if ( class(x) != "data.frame" ) {
                stop("'x' must be a data frame", class.= FALSE)
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_gg_ens = useMart("ensembl", dataset="ggallus_gene_ensembl")
        mart_hsa_ens = useMart("ensembl", dataset="hsapiens_gene_ensembl")
        # Create ortholog table
        ortho_df <- getLDS(attributes=c("ensembl_gene_id","hsapiens_homolog_orthology_confidence"),
                           filters="ensembl_gene_id", values = rownames(x), mart=mart_gg_ens,
                           attributesL=c("hgnc_symbol"), martL=mart_hsa_ens) # Use biomart to get orthologs
        # Filter out any low confidence orthologs and any genes that are not one-to-one orthologs in both directions
        ortho_df <- ortho_df[ortho_df$Human.orthology.confidence..0.low..1.high. == '1',]
        ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
        ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
        
        # Assumes that 'x' has ensembl chicken gene ID's as rownames
        # Ensure that only chicken genes appear in the matrix
        x <- x[startsWith(rownames(x), "ENSGALG"),]
        # Assign the HUGO symbols to a new column
        x$hugo <- NA
        x$chick_ens <- rownames(x)
        
        # Collect the hugo symbols
        ###################################
        get_hugo <- function(x) {
                for ( chick_ens_id in x$chick_ens ) {
                        # Collect the hugo symbol
                        hugo_sym <- ortho_df[ortho_df$Gene.stable.ID == chick_ens_id,]$HGNC.symbol
                        # Assign the hugo symbol to the appropriate column in 'x'
                        x <- mutate(x, hugo = ifelse(x$chick_ens == chick_ens_id, hugo_sym, hugo))
                }
                x
        }
        ###################################
        # Collect the hugo symbols
        x <- get_hugo(x)
        # Remove rows with NA values in hugo field
        x <- x[!is.na(x$hugo),]
        # Set rownames to hugo symbols
        rownames(x) <- x$hugo
        # Remove extra annotation
        drops <- c("hugo", "chick_ens")
        x <- x[,!(names(x) %in% drops)]
        x
}
#########################################################################################

# Capture the Date
date <- format.Date( Sys.Date(), '%Y%m%d' )

#' ## Load Data
#' 
#' Data consists of gene count matrix of chicken genes (only those orthologous 
#' to human included) across normal and tumor samples. Samples are also annotated
#' with relevent phenotypic data. Files have been formatted for easy pipeline comparison via 
#' [BioJupies](https://biojupies.cloud/). Note: load 20190823_biojupiesmatrix_steep.txt into [BioJupies](https://biojupies.cloud/) if you'd like to use gene symbols.
#' 
#' ##### Count Matrix Prefiltering:
#' 
#' ##### Pheno data file annotation:
#' * `sample`: Sample ID
#' * `tissue`: Normal vs. Tumor
#'     * `NC`: Normal Control
#'     * `Tu`: Tumor
#' * `sex`: The sex of the bird
#'     * `M`: Male
#'     * `F`: Female
#' * `ikaros`: Status of IKZF1 mutation or expression (mutation takes presidence)
#'     * `WT`: Wildtype
#'     * `MUT`: Nonsynonymous (always in-frame in this cancer) mutation in IKZF1 coding region
#'     * `LOW`: Comparitevely low IKZF1 expression (compared to normal T cell counts)
#' * `ikaros_func`: Whether Ikaros activity is inferred to remain functional.
#'     * `FUNC`: Function maintained (inferred)
#'     * `LOF`: Loss of DNA-binding domain and tumor suppressor function (inferred)
#' * `mut`: IKZF1 Mutation status (similar to "ikaros" column)
#'     * `NO_MUT`: No nonsynonymous mutation detected in IKZF1
#'     * `MUT`: Nonsynonymous mutation detected in IKZF1
#' * `replicate`: Whether sample is a biological replicate (`YES`/`NO`)
#' * `used`: Adjusted purity estimate from linear model
#' * `purity`: Purity estimate relative to heterozygous (presumably truncal) IKZF1 variant allele frequency

#+ Load the Data

# Load the data
# Pheno data
status <- read.table(file = './data/20190822_metainfo_steep.txt', header = TRUE,
                     row.names = NULL, sep = '\t')

# Raw Counts
all.data <- read.csv(file = './data/20190823_rawcounts_steep.txt', header = TRUE,
                     sep = '\t', check.names=FALSE)

#' #### The Gene Counts (dimhead)
dim(all.data)
head(all.data,3)
#' #### The Pheno Data (dimhead)
dim(status)
head(status, 32)

# Collect female samples
females <- as.character(filter(status, sex == 'F')$sam)
# Collect male samples
males <- as.character(filter(status, sex == 'M')$sam)
# Collect control samples
controls <- as.character(filter(status, tissue == 'NC')$sam)
# Collect samples
samples <- status$sample
normal_samples <- c('017733','017748','017820','017824','017936','017939','017945','017947')
tumor_samples <- as.character(samples[samples %!in% normal_samples])

# "It is absolutely critical that the columns of the count matrix and the rows of the 
# column data (information about samples) are in the same order. DESeq2 will not make 
# guesses as to which column of the count matrix belongs to which row of the column data, 
# these must be provided to DESeq2 already in consistent order." ~ Mike Love
all.data.m <- as.matrix(all.data)
rownames(status) <- status$sample 
all(rownames(status) %in% colnames(all.data.m))
all(rownames(status) == colnames(all.data.m))
all.data.m <- all.data.m[, rownames(status)]
all(rownames(status) == colnames(all.data.m))

##########################################################################
############ Annotate Data for BioJupies #################################
##########################################################################
# Annotate Data for BioJupies
# Annotate normalized counts
biojupes <- all.data
biojupes$ensembl <- rownames(biojupes)
biojupes$symbol <- mapIds(org.Gg.eg.db, biojupes$ensembl, "SYMBOL", "ENSEMBL")
biojupes <- biojupes[!duplicated(biojupes$symbol),]
biojupes <- biojupes[!is.na(biojupes$symbol),]
rownames(biojupes) <- biojupes$symbol
biojupes <- biojupes %>% dplyr::select(symbol, everything())
biojupes <- biojupes %>% dplyr::select(-ensembl)
write.table(biojupes, file = paste0('./data/',date,'_biojupiesmatrix_steep.txt'), quote = FALSE, row.names = FALSE, sep = '\t')
##########################################################################


#' ## PCA Analysis
#' A principal component analysis before surrogate variable and power adjustments reveals The variance between samples of different sexes and estimated tumor purity was reduced, while the variance between the primary dichotomous relationship of interest—tumors versus normal samples—was increased. We believe that the primary, secondary, and tertiary factors influencing these datasets are the relationship between tumor and normal tissues, tumor purity, and the contamination of gonadal tissue in and around tumors from different gonadal phenotypes (determined by sex), respectively.
#' 

#+ PCA Analysis (Unsupervised model)

# Perform unsupervised clustering
# Build model
count_data <- all.data.m
col_data <- status # To look at all samples 
design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
dds
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]
dds

dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
rld <- vst(dds) #vst and rlog comparable with all samples
#rld <- rlog(dds, blind=FALSE)

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts <- as.data.frame(norm_counts)

# Annotate normalized counts
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Gg.eg.db, counts$ensembl, "SYMBOL", "ENSEMBL")
counts$entrez <- mapIds(org.Gg.eg.db, counts$ensembl, "ENTREZID", "ENSEMBL")
counts$genename <- mapIds(org.Gg.eg.db, counts$ensembl, "GENENAME", "ENSEMBL")
counts$go <- mapIds(org.Gg.eg.db, counts$ensembl, "GO", "ENSEMBL")
counts$path <- mapIds(org.Gg.eg.db, counts$ensembl, "PATH", "ENSEMBL")


#' ### Primary variable of interest (tumors vs normals), unsupervised
#' Let's examine our primary dichotomous relationship of interest--tumors vs normals. How well does this relationship explain the variance in our data? Are there other factors influencing our data?


# PCA Analysis (tissue)
DESeq2::plotPCA(rld, intgroup ="tissue") +
        guides(color=guide_legend(title="Tumor vs Normal"))

#' Our normals cluster quite well but there is quite a lot of variance amongst tumor samples.


#' ### Batch effect and technical artifact investigation
#' All sequencing was performed at the same time, place, and machine.
#' Let's look at the lanes ... do we see lane bias? 

#+ Batch effect and technical artifact investigation, echo=TRUE
# PCA Analysis (lane)
DESeq2::plotPCA(rld, intgroup ="lane") +
        guides(color=guide_legend(title="Lane(s)"))

#' Our design was not ideal, we sequenced half of our controls in 1 lane (lane 7). This could be a potential lane bias to watch out for, but I think it's unlikely given distribution of PCs and how tightly the controls cluster (at least compared to tumors).
#'
#' ### Biological Influence (phenotype data) investigation
#' Let's examine the sex of the birds. These tumors seeded in gonads--perhaps the organs that differ most between sexes. Note: Normal samples consist of normal CD4+ T cells and tumor samples contain neoplastic CD4+ T cells seeded in gonadal tissue.

# PCA Analysis (Sex)
DESeq2::plotPCA(rld, intgroup ="sex") +
        guides(color=guide_legend(title="Sex of Bird"))

#' It seems that sex might be driving a significant amount of variance in the gene expression data.

#' Let's examine the estimated purity of tumors. To estimate tumor purity, we used feature selection (LASSO) and a simple linear model between 2 datasets (DNASeq and RNASeq). MD tumors are usually monoclonal (as were most in this experiment) so we were able to identify truncal mutations (VAFs) and compare those to expression counts of genes with highly conserved expression (e.g. CD4 or IKZF1). Tumor purity analysis not included in this pipeline.

# PCA Analysis (Estimated Tumor Purity)
DESeq2::plotPCA(rld, intgroup ="used") +
        guides(color=guide_legend(title="Estimated Neoplastic Cell Ratio"))

#' Tumor purity status might also be driving a significant amount of variance amongst tumors in the gene expression data.
#'
#' Ikaros mutations were found across tumors of different phenotypes. They do, however, appear in higher purity female tumors.

# PCA Analysis (Ikaros Mutations Status)
DESeq2::plotPCA(rld, intgroup ="mut") +
        guides(color=guide_legend(title="Ikaros Mutations Status"))


#' ## Surrogate Variable Analysis
#' To investigate spurious signals, surrogate variable analysis (SVA) was performed incorporating the primary variable of interest: tumor versus normal. There was not significant evidence of influence from technical factors (i.e. date, lane, and location). 
#' The top 2 factors correlated with: the top surrogate variable were sex (p = 2.991e-07) and tumor purity (p = 0.009172). 


# Perform surrogate variable analysis to show a correlation between hypothetical surrogate variables and purity
##################################################################################################################
############################################################################
############### SVA ########################################################
############################################################################

# Compare tumor vs normal to an unsupervised model. Do any other varibales hide behind the comparison of interest (tumor vs normal?)
## Perform surrogate variable analysis (sva) to look for hidden variation in the data
design = ~ tissue
n.sv <- 2 # Number of surrogate variables
# Perform SVA analysis

dat <- counts(dds, normalized=FALSE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(design, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=n.sv)

# Combine the data into a plot
plot_df <- cbind(colData(dds), as.data.frame(svseq$sv))
plot_df <- as.data.frame(plot_df)
# Exmaine 'purity' relationship
ggplot(plot_df, aes(x=purity, y=V1, color = sex)) +
        geom_point() +
        stat_smooth(method = 'lm')
# Perform linear model and see of there's a correlation
mod <- lm(V1 ~ purity, data=plot_df)
summary(mod)
mod <- lm(V1 ~ sex, data=plot_df)
summary(mod)

# Show that sex was successfully corrected for
# Exmaine the sex relationship
ggplot(plot_df, aes(x=sex, y=V1)) +
        geom_jitter(width=0.1) +
        geom_boxplot(alpha = 0.4)
mod <- aov(V1~sex,data = plot_df)
summary(mod)

#' 
#' ### Adjust Model for Tumor Purity
#' Let's make a few judgement calls: Some of our tumor samples contain very little neoplastic cells--all datasets eventually supported these results (cytogenetics, DNASeq, microarrays, and even photos of tumors prior to collection). The bad news, we wasted a fair amount of sequencing money on 5 tumor samples. The good news, we learned how to avoid this mistake for the future and built models to predict which tumors are most valuable for sequencing based on phenotypic data. 
#' 
#' We will (arbitrarily) remove tumor samples with an estimated tumor purity less than 25%. Sex seemed to drive much of the variance amongst tumors, however, we cannot reasonably attempt to correct for this. `sex` really represents infiltrating gondal tissue in tumors; normal samples were T cells. We may be able to adjust for purity. Let's rerun the PCA when adjusting for purity.

##############################################################################
######## Model Adjustment, PCA, Differential Gene Expression Analysis  #######
##############################################################################

# Pheno data
status <- read.table(file = './data/20190822_metainfo_steep.txt', header = TRUE,
                     row.names = NULL, sep = '\t')
# Raw Counts
all.data <- read.csv(file = './data/20190823_rawcounts_steep.txt', header = TRUE,
                     sep = '\t', check.names=FALSE)

# Adjust the status dataframe to only include samples with purity > 25%
status <- filter(status, purity > 0.25) # TCGA maintains higher purity standards
#status <- filter(status, tissue == 'Tu')

# Collect female samples
females <- as.character(filter(status, sex == 'F')$sam)
# Collect male samples
males <- as.character(filter(status, sex == 'M')$sam)
# Collect control samples
controls <- as.character(filter(status, tissue == 'NC')$sam)
# Collect samples
samples <- status$sam

# Make sure rownames of status equal samples
rownames(status) <- status$sam

# Account for estimated purity in design
DE_design = ~ purity

# Save the formula in string form
form <- as.character(DE_design)

# Perform DESeq analysi with males and females seperately (sex of bird is confounding)
all.data.m <- as.matrix(all.data)
all(rownames(status) %in% colnames(all.data.m))
all.data.s <- all.data.m[, rownames(status)]
all(rownames(status) == colnames(all.data.s))
##### Apply differential gene expression analysis with DESeq2
lfc=0.5 
qval=0.05

# Apply differential gene expression analysis with DESeq2
count_data <- all.data.s
col_data <- status
design = DE_design # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)
# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from: 
#citation("PROPER")
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]
dds
# Ensure the factor levels are properly set
dds$tissue <- relevel(dds$tissue, ref = "NC")

# Perform differential expression analysis
dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
rld <- DESeq2::vst(dds, blind=FALSE)
counts <- as.data.frame(assay(rld))
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Gg.eg.db, rownames(counts), "SYMBOL", "ENSEMBL")
dds <- DESeq(dds)
res <- results(dds)
res
res.thr <- results(dds, lfcThreshold=lfc, alpha=qval)
summary(res.thr)
res.thr$ensembl <- rownames(res.thr)
res.thr$symbol <- mapIds(org.Gg.eg.db, rownames(res.thr), "SYMBOL", "ENSEMBL")
res.thr$entrez <- mapIds(org.Gg.eg.db, rownames(res.thr), "ENTREZID", "ENSEMBL")
res.thr$genename <- mapIds(org.Gg.eg.db, rownames(res.thr), "GENENAME", "ENSEMBL")
res.thr$go <- mapIds(org.Gg.eg.db, rownames(res.thr), "GO", "ENSEMBL")
res.thr$path <- mapIds(org.Gg.eg.db, rownames(res.thr), "PATH", "ENSEMBL")
resSort <- res.thr[order(res.thr$padj),]
num <- as.numeric(table(resSort$padj < qval)[2])
topgenes <- head(rownames(resSort),num)

# PCA Analysis (Cohort)
DESeq2::plotPCA(rld, intgroup ="sex") +
        guides(color=guide_legend(title="Cohort"))

#' The variance between samples of different sexes and estimated tumor purity reduced, while the variance between the primary dichotomous relationship of interest—MD tumors versus CD4+ normal samples increased.
#' 
#' 
#' ### Visulaize Differentially Expressed Genese (Tumors vs Normals)

# Extract significant results
resSig <- head(resSort, n = num)
mat <- assay(rld)[topgenes,]
# Expression is regularized log, then subtracted by rowmeans
mat <- mat - rowMeans(mat)

#remove NA values
mat <- mat[!is.na(rownames(mat)),]

# Extract the annotation information
hm_title <- c("sex","purity",'mut')
df <- as.data.frame(colData(dds)[,hm_title])
df$mut <- as.character(df$mut)
df[df$mut == "NO_MUT",]$mut <- 'Wt'
df[df$mut == "MUT",]$mut <- 'Mutant'
df$mut <- as.factor(df$mut)
names(df) <- c('Sample status','Estimated purity','IKZF1 status')

# Subset samples to make plotting easy
#n500 <- sample(nrow(counts), 500)
#counts_hm <- counts[n500,]
# Remove controls
colnames(mat)
#mat <- mat[,c(9:26)]
dim(mat)

#x <- as.data.frame(mat)
#write.table(x, file = './data/mat_tumors_de_genes_4heatmap.txt', quote = FALSE,
#            row.names = FALSE, sep = '\t')

#############################################################################
########## Quickly generate tumor heatmap of gene expression ################
#############################################################################
# Load the file
# x <- read.csv(file = './data/mat_tumors_de_genes_4heatmap.txt', sep = '\t')
#x <- as.matrix(x)
x <- mat

hm <- pheatmap(x, annotation_col=df, fontsize = 16,
               show_rownames=FALSE,show_colnames=FALSE)
hm
#############################################################################
#############################################################################
#############################################################################


#' #### This was a modest attempt to understand what was driving the variance in our gene expression datasets.



