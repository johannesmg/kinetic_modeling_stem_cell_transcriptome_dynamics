#library(tidyverse)
library(timecourse)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)


args = commandArgs(trailingOnly=TRUE)
eset.file=args[1]
out.file=args[2]


batch.compensation.linear.interpolation <- function(x){
                                        #    xx <- function(y){approx(x=x$t,y=y$expr.bt,xout=y)}
    if (length(x$t)>1){
        xx <- approxfun(x=x$t,y=x$expr.bt)
    }else{
        xx <- approxfun(x=0:144,y=rep(x$expr.bt,length(0:144)))
    }
    xx
}

dummy.fn <- function(t,expr.bt,conc,batch){
    oo <- rep(0,length(batch))
    for (i in unique(batch)){
        xx=t[batch==i & conc==0]
        yy=expr.bt[batch==i & conc==0]
        if (length(xx)>1){
            oo[batch==i] <- approx(x=xx,y=yy,xout=t[batch==i])$y
        }else{
            oo[batch==i] <- approx(x=c(min(t),max(t)),y=rep(yy,2),xout=t[batch==i])$y
        }
    }
    oo
}


error.predictor <- function(x){
    predict(ll,newdata=x)
}


compute.error <- function(value){
    (1/2*sd(value)+1/2*error.predictor(mean(value)))/sqrt(length(value))
}


eset <- read_csv(eset.file,guess_max=2e6)


filter(eset,(!is.na(t) & t>=0),compound =="untreated" | is.na(compound),protocol=="UKN1",grepl("UKN1time",sample)) %>%
    filter(study=="timecourse") %>%
    group_by(gene,t) %>%
    summarise(mn=mean(value),sd=sd(value)) %>%
    dplyr::select(mn:sd) %>%
    ungroup() %>%
    mutate(mn.int=cut(mn,15)) %>%
    group_by(mn.int) %>%
    summarise(mn=mean(mn),sd=mean(sd)) %>%
    as.data.frame() ->  mn.sd.frame

mn.sd.frame[mn.sd.frame$mn<6,"sd"] <- max(mn.sd.frame["sd"])
ll <- loess(sd~mn,data=mn.sd.frame,span=0.75,control=loess.control(surface="direct"))


filter(eset,study=="VPA sat conc",t==144,compound!="untreated"|(is.na(treatment_start)&is.na(treatment_end))) %>%
    # biom_A_untreated is outlier
    filter(!grepl("biomAuntrD11",sample))%>%
    filter(compound %in% c("untreated","VPA")) %>%
    filter(is.na(treatment_start)|treatment_start==0,is.na(treatment_end)|treatment_end==144) %>%
    mutate(conc=concentration) %>%
    mutate(conc=ifelse(is.na(conc),0,round(conc/1e-3,3))) %>%
    mutate(t=144,batch=3) %>%
    ungroup %>%
    dplyr::select(one_of(c("gene","t","conc","value","batch","replicate","sample")))%>%
    ungroup ->
    data.batch.3

data.batch.3 %>%
    unite("t_conc",t,conc) %>%
    group_by(t_conc) %>%
    count


#filter(ungroup(eset.time),as.character(X1) %in% all.genes,first.study==TRUE) %>%
filter(eset,(!is.na(t) & t>=0),compound %in% c("untreated","VPA") | is.na(compound),protocol=="UKN1") %>%
    filter(study=="timecourse") %>%
    mutate(first.study=replicate %in% c("A","B","C","D")) %>%
    filter(first.study==TRUE) %>%
    mutate(conc=0,batch=1) %>%
    dplyr::select(one_of(c("gene","t","conc","value","batch","replicate","sample")))->
    data.batch.1
    

data.batch.1 %>%
    unite("t_conc",t,conc) %>%
    group_by(t_conc) %>%
    count


filter(eset,grepl("Dreser",sample)) %>%
    filter(is.na(treatment_start) |treatment_start==0 & treatment_end==t) %>%
    mutate(conc=ifelse(is.na(concentration),0,concentration/1e-3),batch=2) %>%
    dplyr::select(one_of(c("gene","t","conc","value","batch","replicate","sample")))->
    data.batch.2

data.batch.2 %>%
    unite("t_conc",t,conc) %>%
    group_by(t_conc) %>%
    count





data.all.batches <- rbind(data.batch.1,data.batch.2,data.batch.3)




data.all.batches %>%
    group_by(gene,batch,t,conc) %>%
    summarise(mny=mean(value),sdy=compute.error(value)) %>%
    group_by(gene,t) %>%
    mutate(expr.bt=mny-mny[batch==1]) %>%
    filter(batch!=1) %>%
    group_by(gene) %>%
    mutate(mny=mny-dummy.fn(t,expr.bt,conc,batch)) %>%
    dplyr::select(one_of(c("gene","batch","t","conc","mny","sdy"))) %>%
    ungroup ->
    batch.corrected.data


group_by(batch.corrected.data,gene) %>%
    mutate(mny=ifelse(t==144 &conc>0&batch==3,mny-mny[t==144 & conc==0&batch==3],mny)) %>%
    filter(conc>0) ->
    batch.corrected.data


data.all.batches %>%
    filter(batch==1) %>%
    group_by(gene,batch,t,conc) %>%
    summarise(mny=mean(value),sdy=compute.error(value)) ->
    first.batch.data


all.data.batch.corrected <- bind_rows(first.batch.data,batch.corrected.data) %>% ungroup

all.data.batch.corrected <- filter(all.data.batch.corrected,t!=3)




data.batch.1 %>%
    ungroup %>%
    filter(t!=3) %>%
    unite("t_rep",t,replicate) %>%
    dplyr::select(gene,t_rep,value) %>%
    spread(t_rep,value) ->
    dd0

time.vector <- unique(data.batch.1$t)
time.vector <- time.vector[time.vector!=3]
col.vec <- paste(rep(as.character(sort(time.vector)),4),unlist(lapply(c("A","B","C","D"),rep,length(time.vector))),sep="_")
dd0 <- dd0[,c("gene",col.vec)]
out1 <- mb.long(as.matrix(dd0[,-1]), times=12, reps=rep(4,nrow(dd0)))
names(out1$HotellingT2) <- as.character(dd0$gene)

## setup.mart <- function(){
##   ensMart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org")
## }

## get.biotype.human <- function(genes){
##   ensMart=setup.mart()
##   ensMart=useDataset("hsapiens_gene_ensembl",mart=ensMart)
##   df.0 <- getBM(attributes = c("ensembl_gene_id","gene_biotype"),filters = "ensembl_gene_id", values = genes,mart = ensMart)
##   df.0 <- df.0[!duplicated(df.0$"ensembl_gene_id"),]
##   rownames(df.0) <- df.0[,"ensembl_gene_id"]
##   df.0 <- na.omit(df.0)
##   df.0
## }


## biotype <- get.biotype.human(unique(names(sort(out1$HotellingT2,decreasing=TRUE)[1:3000])))
## biotype <- biotype[biotype$gene_biotype=="protein_coding",]

# TIMESTAMP JUNE 29 2016 // ENSEMBL 84 (MARCH 2016)
load(file="./input_data/biotype_human.RData")

sorted.regulated.genes <- names(sort(out1$HotellingT2,decreasing=TRUE))
sorted.regulated.genes <- sorted.regulated.genes[sorted.regulated.genes %in% biotype[,"ensembl_gene_id"]]

sorted.regulated.genes <- sorted.regulated.genes[1:2000]

all.data.batch.corrected <- filter(all.data.batch.corrected,gene %in% sorted.regulated.genes)

write.csv(all.data.batch.corrected,file=out.file,row.names=FALSE) 

write.table(sorted.regulated.genes,file="./output_data/all_changing_genes.csv",row.names=FALSE,col.names=FALSE)










filter(eset,grepl("Dreser",sample)) %>%
    #0_72 E is outlier
    filter(sample!="Dreser111HGU133Plus2") %>%
    mutate(start=treatment_start,end=treatment_end) %>%
    mutate(conc=ifelse(is.na(concentration),0,concentration/1e-3),batch=2) %>%
    dplyr::select(one_of(c("gene","t","start","end","conc","value","batch","replicate","sample")))->
    data.batch.window

data.batch.window %>%
    unite("t_start_end_conc",t,start,end,conc) %>%
    group_by(t_start_end_conc) %>%
    count

data.batch.window %>%
    group_by(gene,batch,t,start,end,conc) %>%
    summarise(mny=mean(value),sdy=compute.error(value)) ->
    data.batch.window

write.csv(data.batch.window,file="./output_data/vpa_window_treatment_data.csv")





filter(eset,grepl("VPA sat acute",study),t==144,is.na(compound)|compound=="untreated") %$%
    unique(sample) ->
    vpa.acute.study.untreated

filter(eset,grepl("VPA sat acute",study),treatment_start==96,treatment_end==144,compound=="VPA",concentration==6e-04) %$%
    unique(sample) ->
    vpa.acute.study.vpa

filter(eset,sample %in% c(vpa.acute.study.untreated,vpa.acute.study.vpa)) %>%
    mutate(start=treatment_start,end=treatment_end) %>%
    mutate(conc=ifelse(is.na(concentration),0,concentration/1e-3),batch=2) %>%
    dplyr::select(one_of(c("gene","t","start","end","conc","value","batch","replicate","sample"))) %>%
    group_by(gene,batch,t,start,end,conc) %>%
    summarise(mny=mean(value),sdy=compute.error(value)) ->
    data.batch.window.vpa.acute

write.csv(data.batch.window.vpa.acute,file="./output_data/vpa_acute_study_window_treatment_data.csv")



