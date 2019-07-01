#rm(list=ls())
#-- define the pattern of files to be analysed, the file shold end as _HtSeqCount.txt
setwd("/Users/Pooja/Documents/Data-Analysis/Others/FangFang/2018/ribosome_profiling")
count_file="/Users/Pooja/Documents/Data-Analysis/Others/FangFang/2018/ribosome_profiling/count_matrix_bed_coverage.txt"
 
#-- provide metadata file containg file names and their the replicate information as below
 metadata_file <- "metadata.txt"
 
#' Title: deseq_from_prepDE
#' deseq_from_prepDE(count_file,metadata_file, "outfile") 
#' @param dir 
#' @param pattern 
#' @param metadata_file 
#'treatment_set1	treated
# treatment_set2	treated
# control_set1	untreated
# control_set2	untreated
#' @param outfile 
#' @author Pooja Sethiya (yb57662@umac.mo)
#' University of Macau
#' @return
#' @export
#'
#' @examples


eseq_from_prepDE <- function(count_file,metadata_file, outfile){
          
          #--- Load package
          library(tidyverse)
          library(data.table)
          library(DESeq2)
          library(purrr)
          library(GGally)
          library(xlsx)
          
          #--- get count matrix (prepDE output)         
          df_count_matrix <- read_delim(count_file, delim=",", col_names = T)
          condition <- read_delim(metadata_file, delim="\t", col_names = FALSE)
          
          
          #--- filter genes if read count is less than 2 in all columns
     
          df_count_matrix_filtered <- df_count_matrix %>% filter_if(is.numeric, any_vars(. >= 2)) %>% as.tibble(rownames="gene_id")
          
          
          countdata <- as.data.frame(df_count_matrix_filtered) %>% 
                                       column_to_rownames("gene_id")%>%
                                       dplyr::select(condition$X1) %>% # keep sample columns using sampleinfo
                                       as.matrix()
          
          dim(countdata)
          
          colData <- as.data.frame(colnames(countdata))
          colData$condition <- condition$X2[match(colData[,1], condition$X1)]
          
          colnames(colData) <- c("colData", "condition")
          print(ncol(countdata) == nrow(colData))
          
          colData$condition=factor(colData$condition, levels=c("untreated", "treated"))
          print(colData)
          
          #--- make DESeq object
          
          ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = colData,design = ~ condition)
          
          
          
      
          vstcounts <- vst(ddsObj, blind=TRUE)
          print(plotPCA(vstcounts))
          
          dds <- estimateSizeFactors(ddsObj) # run this or DESeq() first
          
          dds$condition <- factor(dds$condition)
       
          #--- Compute DESeq2
          dds  <- DESeq(dds) # for replicates
          
          res <- results(dds)
          cat(summary(res))
          
          #--- get normalised counts
          resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE) 
          
          print(head(resdata))
          
          expression_matrix <- resdata[,c(1,8:ncol(resdata))]
          print(head(expression_matrix))
          
          deseq_matrix <- expression_matrix %>% column_to_rownames("Row.names")
          #plot_correlation
          gg <-  
                    ggpairs(log2(deseq_matrix+1), upper = list(continuous = wrap("cor", size = 12, color="red"))) +
                    theme_bw() +
                    theme(legend.text = element_text(size = 25),
                          legend.title = element_text(size = 20),         
                          axis.title.x = element_text(size = 15),        
                          axis.title.y = element_text(size = 15),
                          axis.text.x = element_text(size = 12),
                          axis.text.y = element_text(size = 12),
                          strip.text = element_text(size = 12, color="black")) 
          
          png(paste(outfile, "_correlation1_plot.png", sep=""),width = 1200, height = 1200, units = "px", pointsize = 12)
          print(gg )
          dev.off()
          
          up_deg <- subset(resdata, resdata$log2FoldChange> 0.6)
          down_deg <- subset(resdata, resdata$log2FoldChange< -0.6)
       
          write.xlsx(resdata,sheetName = "deseq_output", file=paste(outfile, "_deseq_output.xlsx", sep=""),row.names = F, append=F)
          write.xlsx(up_deg,sheetName = "up_DEG",file=paste(outfile, "_deseq_output.xlsx", sep=""),row.names = F, append=T)
          write.xlsx(down_deg, sheetName ="down_DEG",file=paste(outfile, "_deseq_output.xlsx", sep=""),row.names = F, append=T)
          
}













         