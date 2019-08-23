#'---
#' title: "RNASeq Workflow: Specific Usage"
#' author: Alec Steep
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' # Bulk RNA Seq Analysis of Marek's disease lymphomas seeded in gonadal tissue
#' 
#' ## Rudimentary Differential Gene Expression Workflow
#' * Explanation of Samples and Data
#' * Download Data:
#'     * Gene count matrix & metadata file
#'         * Files formatted for pipeline comparison with [BioJupies](https://biojupies.cloud/)
#' * 
#' 
#' ## Explanation of Samples and Data 
#' ### Explaination of Biological Samples and Sequencing
#' Biological Samples:
#' Marek's disease (MD) lymphomas seeded in the gonads of males and females 
#' of F1 birds--cross of parental lines 6 (resistent to MD) and 7 (suspectible to 
#' MD)--infected at hatch with 1000 pfu of Marek's disease virus (MDV) (strain JM102W). 
#' Birds were euthanized and tumors collected when birds were moribund or if they 
#' reached eight weeks of age (birds challenged at hatch). RNA was extracted via 
#' the miTNeasy mini kit (Qiagen). Sequencing was performed in 1 batch; 26 gonadal 
#' tumor samples from 22 tumors from 22 birds and 8 samples of isolated CD4+ spleenic 
#' T-cells from 8 uninfected birds on 05/28/2015. Experimental design decision:
#' normal CD4+ T cell populations cannot be obtained from MD-infected birds;
#' therefore, normal CD4+ cell populations were obtained from uninfected birds
#' of 6x7 F1 cross. RNA-Seq libraries were prepared using the NuGen Ovation 
#' Single Cell RNA-Seq System (stranded). All samples underwent sequencing to 
#' produce 125bp paired-end reads via Illumina HiSeq machines at the Michigan 
#' State University RTSF Genomics Core. Project designed and performed by 
#' Hans Cheng (Hans.Cheng@ars.usda.gov) and Alec Steep (alec.steep@gmail.com).
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
#' ## Download Data
#' 
#' Data consists of gene count matrix of chicken genes (only those orthologous 
#' to human included) across normal and tumor samples. Samples are also annotated
#' with relevent metadata. Files have been formatted for pipeline comparison with 
#' [BioJupies](https://biojupies.cloud/).
#' 
#' ##### Count Matrix Prefiltering:
#' 
#' ##### Metadata file annotation:
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

#+ Setup Environment,

############################################################
##### Dependencies #########################################
############################################################

# Load dependencies
#################
#BiocManager::install("SeqGSEA", version = "3.8")
pacs...man <- c("GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2",
                "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr",
                "org.Gg.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt","feather",
                "PROPER","SeqGSEA",'purrr','BioInstaller')

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
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}

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

#+ Download Data

# Pull data from git repository (NTD)





#' ## PCA Analysis
#' A principal component analysis is evident in Figure “PCA” before surrogate variable and power adjustments.  The variance between samples of different sexes and estimated tumor purity was reduced, while the variance between the primary dichotomous relationship of interest—tumors versus normal samples—was increased. We believe that the primary, secondary, and tertiary factors influencing these datasets are the relationship between tumor and normal tissues, tumor purity, and the contamination of gonadal tissue in and around tumors from different gonadal phenotypes (determined by sex), respectively.
#' 

#+ PCA Analysis


#' 
#' 
#' 
#' Normalization for estimated tumor purity was performed and samples with estimated purity less than 0.57 were removed. The variance between samples of different sexes and estimated tumor purity reduced, while the variance between the primary dichotomous relationship of interest—MD tumors versus CD4+ normal samples increased.
#' 
#' 
#' 

#' ## Surrogate Variable Analysis
#' To investigate spurious signals, surrogate variable analysis (SVA) was performed incorporating the primary variable of interest: tumor versus normal. There was not significant evidence of influence from technical factors (i.e. date, lane, and location). The top 2 factors correlated with one surrogate variable were sex (p = 0.0000000349) and tumor purity (p = 0.009465). IKZF1 mutation status demonstrated association with the second surrogate variable, suggesting the possibility of tumor subtypes based on IKZF1 mutation status (p-value: 0.0303). 
#' 
#' 
#' 
#' 
#' 




#' ## Generate document body from comments
#' All the features from markdown and markdown supported within .Rmd documents, I was able to
#' get from within R scripts.  Here are some that I tested and use most frequently:  
#' 
#' * Smart comment fomatting in your R script generate the body and headers of the document
#'     * Simply tweak your comments to begin with `#'` instead of just `#`  
#' * Create markdown headers as normal: `#' #` for h1, `#' ##` for h2, etc.
#' * Add two spaces to the end of a comment line to start a new line (just like regular markdown)
#' * Add `toc: true` to YAML frontmatter to create a table of contents with links like the one at the 
#' top of this page that links to h1, h2 & h3's indented like so:
#'     * h1
#'         * h2
#'             * h3
#' * Modify YAML to change syntax highlighting style (I'm using zenburn), author, title, theme, and all the good stuff
#' you're used to setting in Rmd documents.
#' * Sub-bullets like the ones above are created by a `#' *` with 4 spaces per level of indentation.
#' * Surround text with `*` to *italicize*  
#' * Surround text with `**` to **bold**  
#' * Surround text with `***` to ***italicize & bold***  
#' * Skip lines with `#' <br>`
#' * Keep comments in code, but hide from printing in report with `#' <!-- this text will not print in report -->`  
#' * Add hyperlinks:
#'     * [Rmarkdown cheatsheet](http://rmarkdown.rstudio.com/RMarkdownCheatSheet.pdf)
#'     * [Rmarkdown Reference Guide](http://rmarkdown.rstudio.com/RMarkdownReferenceGuide.pdf)
#'     * [Compiling R notebooks from R Scripts](http://rmarkdown.rstudio.com/r_notebook_format.html)
#' 

