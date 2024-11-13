#install required package
BiocManager::install(c("airway","DESeq2","EnhancedVolcano"))
install.packages("ggthemes")
#load required packages
library(tidyverse)
library(airway)
library(DESeq2)
#get the data
data("airway")
airway
#counts data
counts_data<-assay(airway)
#metadata (call data)
col_data<-colData(airway)
view(col_data)
as.data.frame(colData(airway))
#making sure rownames in colData(metadata) matches to the column names in count_data
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data)%in% rownames(col_data))
#are they in same order
all(colnames(counts_data)== rownames(col_data))
#contruct DESeqDataSetFromMatrix data object
dds<-DESeqDataSetFromMatrix(
  countData=counts_data, colData = col_data,
  design = ~dex #condition
  )

#prefiltering
#removing rows with low gene counts/keep rows that at least 10 reads total
rowSums(counts(dds))
keep<-rowSums(counts(dds))>=10
#keep

#reference category
dds$dex<-relevel(dds$dex, ref = "untrt")
#perform differential gene expression analysis
dds<-DESeq(dds)
#save the results
res<-results(dds)
res_df<-as.data.frame(res)

#exploring results
summary(res)
#working with alpha (significance level)
#if alpha=0.1 or 10%~90% CI (default)
#if alpha=0.01 or 1%~99% CI (custom significance level)
#if alpha=0.05 or 1%~95% CI (custom significance level)

res_0.01<-results(dds, alpha = 0.01)
#get summary
summary(res_0.01)

res_0.05<-results(dds, alpha = 0.05)
summary(res_0.05)
#results name
resultsNames(dds)

#contrast
contrast_res<-results(dds,contrast = c("dex","trt","untrt"))
contrast_res
summary(contrast_res)

#to converts data as data frame


#export the data
res_0.01<-results(dds, alpha = 0.01)
res_0.01_df<-as.data.frame(res_0.01)
write.csv(res_0.01_df, "airway.csv", row.names = F)
