library(GEOquery)
library(GSEABase)
library(reshape2)
library(tidyverse)
library(org.Hs.eg.db)


library(tidyverse)
library(magrittr)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
symbol.conversion.file=args[1]
out.file=args[2]



load(symbol.conversion.file)

gset <- getGEO("GSE45223", GSEMatrix =TRUE, getGPL=T)

mica.eset=exprs(gset[[1]])
mica.exprs <- bind_cols(tibble(ID_REF=rownames(mica.eset)),as_tibble(mica.eset))
mica.pheno <- as_tibble(pData(gset[[1]]))

group_by(mica.pheno,source_name_ch1) %>%
    do(ll=unique(.$geo_accession)) %$%
    setNames(ll,source_name_ch1) ->
    condition.list


nn <- mica.exprs$ID_REF
nn.2 <- ga[["symbol.ensg"]][nn]
is.na.or.duplicate <- !is.na(nn.2) & !duplicated(nn.2)
mica.exprs <- mica.exprs[is.na.or.duplicate,]
mica.exprs$gene <- nn.2[is.na.or.duplicate]

dm.fn <- function(x){
    rv <- apply(mica.exprs[,x],1,mean)
    names(rv) <- mica.exprs$gene
    rv
}

mn.sd <- function(x){
    rv.mn <- apply(mica.exprs[,x],1,mean)
    rv.sd <- apply(mica.exprs[,x],1,sd)
    rv <- data.frame(mn=rv.mn,sd=rv.sd)
    rownames(rv) <- mica.exprs$gene
    rv
}

select.fn <- function(cond.1,cond.2,thrsh){
    z.1 <- z.value.dist(cond.1)
    z.2 <- z.value.dist(cond.2)
    mean.1 <- apply(mica.exprs[,condition.list[[cond.1]]],1,mean)
    mean.2 <- apply(mica.exprs[,condition.list[[cond.2]]],1,mean)
    df=mean.1-mean.2
    names(df) <- mica.exprs$gene
    df <- df/mica.error.fn(1/2*(mean.1+mean.2))
    pos <- names(df)[df>thrsh]
    pos.2 <- names(df)[df< -thrsh]
    fdr.1 <- sum(z.1>thrsh)
    fdr.2 <- sum(z.2>thrsh)
    list(pos,pos.2,fdr.1/(length(pos)+length(pos.2)),fdr.2/(length(pos)+length(pos.2)))
}

select.fn.df <- function(cond.1,cond.2,thrsh){
    z.1 <- z.value.dist(cond.1)
    z.2 <- z.value.dist(cond.2)
    mean.1 <- apply(mica.exprs[,condition.list[[cond.1]]],1,mean)
    mean.2 <- apply(mica.exprs[,condition.list[[cond.2]]],1,mean)
    df=mean.1-mean.2
    names(df) <- mica.exprs$gene
    df <- df/mica.error.fn(1/2*(mean.1+mean.2))
    df <- enframe(df)
    df <- arrange(df,-value)
    df$symbol=ga[["ensg.symbol"]][df$name]
    df
}

error.normalized.difference <- function(cond.1,cond.2,thrsh){
    z.1 <- z.value.dist(cond.1)
    z.2 <- z.value.dist(cond.2)
    mean.1 <- apply(mica.exprs[,condition.list[[cond.1]]],1,mean)
    mean.2 <- apply(mica.exprs[,condition.list[[cond.2]]],1,mean)
    df=mean.1-mean.2
    names(df) <- mica.exprs$gene
    df <- df/mica.error.fn(1/2*(mean.1+mean.2))
}

mica.mean <- lapply(condition.list,dm.fn)
#mica.mean <- lapply(mica.mean,function(x){x[names(x) %in% all.genes]})
mica.mean.sd <- mclapply(condition.list,mn.sd,mc.cores=12)

do.call("rbind",mica.mean.sd) %>%
    mutate(mn.int=cut(mn,15)) %>%
    group_by(mn.int) %>%
    summarise(mn=mean(mn),sd=mean(sd)) %>%
    as.data.frame() ->  mica.mean.sd

mica.error.fn <- function(x){
    predict(mica.error.model,x)
}

z.value.dist <- function(x){
    x <- condition.list[[x]]
    rv.mn <- apply(mica.exprs[,x],1,mean)
    rv.sd <- apply(mica.exprs[,x],1,sd)
    zvals <- rv.sd/mica.error.fn(rv.mn)
}

mica.mean.sd[mica.mean.sd$mn<10,"sd"] <- max(mica.mean.sd$sd)
mica.error.model <- loess(sd~mn,data=mica.mean.sd,span=0.75,control=loess.control(surface="direct"))



genes.wnt.up <- select.fn("WA09 embryonic stem cells, NC, day 3","WA09 embryonic stem cells, DSi, day 3",3.5)[[1]]
genes.wnt.down <- select.fn("WA09 embryonic stem cells, DSi, day 3","WA09 embryonic stem cells, NC, day 3",3.5)[[1]]
genes.wnt.smad.up <- select.fn("WA09 embryonic stem cells, NC, day 11","WA09 embryonic stem cells, DSi, day 11",3.5)[[1]]
genes.wnt.smad.down <- select.fn("WA09 embryonic stem cells, DSi, day 11","WA09 embryonic stem cells, NC, day 11",3.5)[[1]]

genes.wnt.frame <- select.fn.df("WA09 embryonic stem cells, NC, day 3","WA09 embryonic stem cells, DSi, day 3",3.5)

error.normalized.difference("WA09 embryonic stem cells, NC, day 3","WA09 embryonic stem cells, DSi, day 3",3.5) %>%
    enframe %>%
    rename(Wnt=value) %>%
    left_join(error.normalized.difference("WA09 embryonic stem cells, NC, day 11","WA09 embryonic stem cells, DSi, day 11",3.5) %>% enframe,by="name") %>%
    rename(SMAD=value) %>%
    rename(gene=name) %>%
    mutate(symbol=ga[["ensg.symbol"]][gene]) ->
    z.value.frame


gene.sets <- list(SMAD_Wnt_up=genes.wnt.smad.up,SMAD_Wnt_dn=genes.wnt.smad.down,Wnt_up=genes.wnt.up,Wnt_dn=genes.wnt.down,Wnt_frame=genes.wnt.frame,z_val_frame=z.value.frame)
save(gene.sets,file=out.file)


z.value.frame %>%
    select(gene,symbol,Wnt) %>%
    filter(abs(Wnt)>3.5) %>%
    arrange(Wnt) %>%
    rename(z_value_wnt_induction=Wnt) %>%
    write_csv("./output_data/supplementary_table_wnt_targets.csv")
