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
no.cores=as.numeric(args[5])
out.file=args[6]


source("gsea_functions.R")

gsea.frame.in <- read_csv(gsea.in.file)

gsea.frame.in %>%
    group_by(sample) %>%
    do(ll=deframe(.[c("gene","value")])) %$%
    setNames(ll,sample) ->
    gsea.list.in

aa <- getGmt(set.definition.file)
gene.set.list <- geneIds(aa)

gsea.out <- mclapply(1:length(gsea.list.in),function(x){print(sprintf("%02f",x/length(gsea.list.in)));rv <- own.gsea(data.vector=gsea.list.in[[x]],ListGSC=list(broad=gene.set.list),min.set.size=minsetsize,nperm=noreps,no.cores.gsea=1);rv},mc.cores=no.cores)
names(gsea.out) <- names(gsea.list.in)

gsea.out.frame <- lapply(gsea.out,gsea.to.frame)
gsea.out.frame <- melt(gsea.out.frame,id.vars=1:4)
gsea.out.frame$min.set.size=minsetsize
gsea.out.frame$no.reps=noreps
gsea.out.frame$in.file=gsea.in.file
gsea.out.frame$set.definition.file=set.definition.file

write_csv(gsea.out.frame,path=paste(out.file))

