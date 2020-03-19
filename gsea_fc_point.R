#source("~/projektentwicklung/sysdt/cluster_compatible_kinetic_model/new_functions.R")

library(dplyr)
library(tidyr)
library(readr)
library(magrittr)

library(tibble)
library(parallel)
library(GSEABase)
library(HTSanalyzeR)
library(reshape2)


args = commandArgs(trailingOnly=TRUE)
gsea.in.file=args[1]
set.definition.file=args[2]
minsetsize=as.numeric(args[3])
noreps=as.numeric(args[4])
no.cores.gsea=as.numeric(args[5])
current.type=args[6]
current.time=args[7]
current.conc=args[8]
out.file=args[9]


source("gsea_functions.R")

gsea.frame.in <- read_csv(gsea.in.file)

current.sample=paste(current.type,current.time,current.conc,sep="_")

gsea.frame.in %>%
    filter(sample==current.sample) ->
    gsea.data.in

gsea.data.in <- deframe(gsea.data.in[c("gene","value")])

aa <- getGmt(set.definition.file)
gene.set.list <- geneIds(aa)

#print(head(gsea.data.in))

gsea.out <- own.gsea(data.vector=gsea.data.in,ListGSC=list(gene_set=gene.set.list),min.set.size=minsetsize,nperm=noreps,no.cores.gsea=no.cores.gsea)


gsea.out.frame <- gsea.to.frame(gsea.out)
gsea.out.frame$sample=current.sample
gsea.out.frame$min.set.size=minsetsize
gsea.out.frame$no.reps=noreps
gsea.out.frame$set.definition.file=set.definition.file

write_csv(gsea.out.frame,path=out.file)

