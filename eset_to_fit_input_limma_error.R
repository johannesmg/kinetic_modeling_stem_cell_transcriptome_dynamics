library(timecourse)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(readr)
library(limma)

args = commandArgs(trailingOnly=TRUE)
eset.file=args[1]
output.file=args[2]

eset.annotated=read_csv(eset.file,guess_max=1e6)

eset.annotated %>%
    filter(t>=0) %>%
    filter(t!=3) %>%
    distinct(study,t,concentration,treatment_start,treatment_end,replicate) %>%
    filter(is.na(treatment_end)|treatment_end==t) %>%
    filter(is.na(treatment_start)|treatment_start==0) %>%
    arrange(study,t) %>%
    group_by(study,t,concentration,treatment_start,treatment_end) %>%
    summarise(nn=n()) %>%
    write_csv(path="./output_data/overview_samples.csv")

eset.annotated %>%
    filter(t>=0) %>%
    filter(t!=3) %>%
    distinct(study,t,concentration,treatment_start,treatment_end,replicate) %>%
    filter(!is.na(treatment_end)) %>%
    filter(!is.na(treatment_start)) %>%
    filter(study=="timecourse") %>%
    group_by(study,t,concentration,treatment_start,treatment_end) %>%
    summarise(nn=n())


eset.annotated %>%
    filter(sample!="Dreser111HGU133Plus2") ->
    eset.annotated

# Endpoint concentration 8*3 + timecourse utr 12*4+ timecourse vpa 2*3+4*2 = 86 samples


eset.annotated %>%
    select(sample,t,replicate,compound,treatment_start,treatment_end,concentration,batch) %>%
    distinct ->
    sample.translation.table



## tibble(sample=colnames(eset.reduced)) %>%
##     left_join(sample.translation.table) %>%
sample.translation.table %>%
    filter(!is.na(batch)) %>%
    unite("condition",t,concentration,treatment_start,treatment_end,batch) %>%
    select(sample,condition) %>%
    mutate(condition=str_replace(condition,"-72","minus72")) %>%
    mutate(foo=1) %>%
    spread(condition,foo) ->
    d0

d0[is.na(d0)]=0

dd0 <- as.data.frame(d0[,-1])
rownames(dd0)=d0$sample


eset.annotated %>%
    select(gene,sample,value) %>%
    spread(sample,value) ->
    eset.spread

eset.reduced=as.matrix(eset.spread[,-1])
rownames(eset.reduced)=eset.spread$gene

sn <- lapply(colnames(dd0),function(x){rownames(dd0)[dd0[,x]==1]})
names(sn) <- colnames(dd0)
fits <- lapply(sn,function(x){lmFit(eset.reduced[,x])})
eb <- lapply(fits,eBayes)

lapply(eb,function(x){head(x$s2.post)})

## tibble(sample=colnames(dd0),m.sd=unlist(lapply(eb,function(x){mean(x$s2.post)}))) %>%
##     arrange(-m.sd) %>%
##     as.data.frame


eb.error.frame=as_tibble(do.call("cbind",lapply(eb,"[[","s2.post")))

bind_cols(tibble(gene=names(eb[[1]]$s2.post)),eb.error.frame) %>%
    gather("sample","value",-gene) ->
    eb.error.frame

no.replicates=as_tibble(do.call("cbind",lapply(sn,length))) %>%
    gather("sample","length")

# standard error of the mean = \sigma/sqrt(n)

eb.error.frame %>%
    left_join(no.replicates) %>%
    mutate(value=sqrt(value)) ->
    eb.error.frame.sd.not.sem


eb.error.frame %>%
    left_join(no.replicates) %>%
    mutate(value=sqrt(value)/sqrt(length)) ->
    eb.error.frame


eset.annotated %>%
    ungroup %>%
    dplyr::filter(t!=3) %>%
# 144_72_144_E is outlier and should apparently be in reality 144_0_96!    
    filter(sample!="Dreser111HGU133Plus2") %>%
    select(gene,t,concentration,treatment_start,treatment_end,replicate,batch,value) %>%
    unite("sample",t,concentration,treatment_start,treatment_end,batch) %>%
    group_by(gene,sample) %>%
    summarise(y=mean(value)) ->
    eset.mean

eset.mean %>%
    left_join(eb.error.frame.sd.not.sem) %>%
    filter(!sample=="-72_NA_NA_NA_1") %>%
#    filter(!is.na(batch)) %>%
    rename(mny=y,sdy=value) %>%
    select(-length) ->
    eset.mean.sd.not.sem


eset.mean %>%
    left_join(eb.error.frame) %>%
    filter(!sample=="-72_NA_NA_NA_1") %>%
#    filter(!is.na(batch)) %>%
    rename(mny=y,sdy=value) %>%
    select(-length) ->
    eset.mean


## eset.mean %>%
##     ungroup %>%
##     separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
##     mutate(t=as.numeric(t)) %>%
##     filter(conc=="NA") %>%
##     mutate(batch=as.factor(batch)) %>%
##     group_by(t,batch) %>%
##     summarise(msdy=mean(sdy)) %>%
##     ggplot(aes(x=t,y=msdy,group=batch,colour=batch))+
##     geom_point()+
##     geom_line()


## eset.mean %>%
##     ungroup %>%
##     separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
##     mutate(t=as.numeric(t)) %>%
##     mutate(conc=as.numeric(conc)) %>%
##     filter(t=="144") %>%
##     filter(batch==3) %>%
##     group_by(t,conc) %>%
##     summarise(msdy=mean(sdy)) %>%
##     ggplot(aes(x=conc,y=msdy,group=1))+
##     geom_point()+
##     geom_line()

## eset.mean %>%
##     ungroup %>%
##     group_by(sample) %>%
##     summarise(msdy=mean(sdy)) %>%
##     ggplot(aes(x=sample,y=msdy,group=1))+
##     geom_point()+
##     theme(axis.text.x=element_text(angle=90))


## eset.mean %>%
##     ungroup %>%
##     separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
##     mutate(t=as.numeric(t)) %>%
##     filter(conc=="NA") %>%
##     mutate(mny=cut(mny,50,include.lowest=T)) %>%
##     group_by(mny) %>%
##     summarise(msdy=mean(sdy)) %>%
##     ggplot(aes(x=mny,y=msdy,group=1))+
##     geom_line()



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

eset.mean %>%
    ungroup %>%
    separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
    mutate(batch=as.numeric(batch)) %>%
    filter(!is.na(batch)) %>%
    mutate(treatment_start=as.numeric(treatment_start)) %>%
    mutate(treatment_end=as.numeric(treatment_end)) %>%
    mutate(conc=as.numeric(conc)) %>%
    mutate(conc=ifelse(is.na(conc),0,conc/1e-3)) %>%
    mutate(t=as.numeric(t)) %>%
    filter(is.na(treatment_start)|treatment_start==0) %>%
    filter(is.na(treatment_end)|treatment_end==t)  %>%
    mutate(conc=ifelse(is.na(conc),0,conc)) %>%
    group_by(gene,t) %>%
    mutate(expr.bt=mny-mny[batch==1]) %>%
#    filter(batch!=1) %>%
    group_by(gene) %>%
    mutate(mny=mny-dummy.fn(t,expr.bt,conc,batch)) %>%
    dplyr::select(one_of(c("gene","batch","t","conc","mny","sdy"))) %>%
    mutate(mny=ifelse(t==144 &conc>0&batch==3,mny-mny[t==144 & conc==0&batch==3],mny)) %>%
#    mutate(sdy=sqrt(sdy)) %>%
    filter(conc>0|batch==1) %>%
    ungroup ->
    batch.corrected.data

eset.mean.sd.not.sem %>%
    ungroup %>%
    separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
    mutate(batch=as.numeric(batch)) %>%
    filter(!is.na(batch)) %>%
    mutate(treatment_start=as.numeric(treatment_start)) %>%
    mutate(treatment_end=as.numeric(treatment_end)) %>%
    mutate(conc=as.numeric(conc)) %>%
    mutate(conc=ifelse(is.na(conc),0,conc/1e-3)) %>%
    mutate(t=as.numeric(t)) %>%
    filter(is.na(treatment_start)|treatment_start==0) %>%
    filter(is.na(treatment_end)|treatment_end==t)  %>%
    mutate(conc=ifelse(is.na(conc),0,conc)) %>%
    group_by(gene,t) %>%
    mutate(expr.bt=mny-mny[batch==1]) %>%
#    filter(batch!=1) %>%
    group_by(gene) %>%
    mutate(mny=mny-dummy.fn(t,expr.bt,conc,batch)) %>%
    dplyr::select(one_of(c("gene","batch","t","conc","mny","sdy"))) %>%
    mutate(mny=ifelse(t==144 &conc>0&batch==3,mny-mny[t==144 & conc==0&batch==3],mny)) %>%
#    mutate(sdy=sqrt(sdy)) %>%
    filter(conc>0|batch==1) %>%
    ungroup ->
    batch.corrected.data.sd.not.sem

## batch.corrected.data %>%
##     filter(gene=="ENSG00000000003") %>%
##     arrange(batch,t,conc) %>%
##     as.data.frame


## eset.mean %>%
##     ungroup %>%
##     separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
##     mutate(batch=as.numeric(batch)) %>%
##     filter(!is.na(batch)) %>%
##     filter(gene=="ENSG00000000003") %>%
##     arrange(batch,t,conc) %>%
##     as.data.frame

## batch.corrected.data %>%
##     filter(gene==ga[["symbol.ensg"]]["DAZL"]) %>%
##     filter(batch%in% c(1,2))%>%
##     filter(!(batch==2 & conc==0)) %>%
##     ggplot(aes(x=t,y=mny,group=batch,colour=batch))+
##     geom_line()+
##     geom_point()

## batch.corrected.data %>%
##     group_by(gene) %>%
##     ## mutate(mny=ifelse(t==144 &conc>0&batch==3,mny-mny[t==144 & conc==0&batch==3],mny)) %>%
##     ## mutate(sdy=sqrt(sdy)) %>%
##     filter(gene==ga[["symbol.ensg"]]["OLFM3"]) %>%
## #    filter(conc>0) %>%
##     filter(batch%in% c(3))%>%
## #    filter(!(batch==2 & conc==0)) %>%
##     ggplot(aes(x=conc,y=mny,ymin=mny-sdy,ymax=mny+sdy,group=batch,colour=batch))+
##     geom_line()+
##     geom_point()+
##     geom_errorbar()

#save.image("ws_fit_input.RData")

eset.annotated %>%
    filter(batch==1) %>%
    select(gene,t,replicate,value) %>%
    ungroup %>%
    filter(t>=0) %>%
    filter(t!=3) %>%
    unite("t_rep",t,replicate) %>%
    dplyr::select(gene,t_rep,value) %>%
    spread(t_rep,value) ->
    dd0

time.vector <- unique(eset.annotated$t)
time.vector <- time.vector[time.vector!=3]
time.vector <- time.vector[time.vector!= -72]
col.vec <- paste(rep(as.character(sort(time.vector)),4),unlist(lapply(c("A","B","C","D"),rep,length(time.vector))),sep="_")
dd0 <- dd0[,c("gene",col.vec)]
out1 <- mb.long(as.matrix(dd0[,-1]), times=12, reps=rep(4,nrow(dd0)))
names(out1$HotellingT2) <- as.character(dd0$gene)


load(file="./input_data/biotype_human.RData")
sorted.regulated.genes <- names(sort(out1$HotellingT2,decreasing=TRUE))
sorted.regulated.genes <- sorted.regulated.genes[sorted.regulated.genes %in% biotype[,"ensembl_gene_id"]]
sorted.regulated.genes <- sorted.regulated.genes[1:2000]

batch.corrected.data <- filter(batch.corrected.data,gene %in% sorted.regulated.genes)
batch.corrected.data.sd.not.sem <- filter(batch.corrected.data.sd.not.sem,gene %in% sorted.regulated.genes)

#old.data=read_csv("~/projektentwicklung/sysdt/vpa_modeling_figures_consolidated/output_data/fit_data_in.csv")

write.csv(batch.corrected.data,file=output.file,row.names=FALSE) 
write.csv(batch.corrected.data.sd.not.sem,file=paste(output.file,"not_sem_corrected.csv",sep=""),row.names=FALSE) 


write.table(sorted.regulated.genes,file="./output_data/all_changing_genes_limma.csv",row.names=FALSE,col.names=FALSE)



eset.mean %>%
    ungroup %>%
    separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
    mutate(batch=as.numeric(batch)) %>%
    filter(batch==2) %>%
    mutate(start=as.numeric(treatment_start)) %>%
    mutate(end=as.numeric(treatment_end)) %>%
    mutate(conc=as.numeric(conc)) %>%
    mutate(conc=ifelse(is.na(conc),0,conc/1e-3)) %>%
    mutate(t=as.numeric(t)) %>%
    dplyr::select(one_of(c("gene","t","start","end","conc","mny","batch"))) ->
    data.batch.window


write.csv(data.batch.window,file="./output_data/vpa_window_treatment_data.csv")


eset.mean %>%
    ungroup %>%
    separate(sample,into=c("t","conc","treatment_start","treatment_end","batch"),sep="_") %>%
    mutate(batch=as.numeric(batch)) %>%
    filter(batch==4) %>%
    mutate(start=as.numeric(treatment_start)) %>%
    mutate(end=as.numeric(treatment_end)) %>%
    mutate(conc=as.numeric(conc)) %>%
    mutate(conc=ifelse(is.na(conc),0,conc/1e-3)) %>%
    filter(conc %in% c(0,0.6)) %>%
    mutate(t=as.numeric(t)) %>%
    dplyr::select(one_of(c("gene","t","start","end","conc","mny","batch"))) ->
    data.batch.window.vpa.acute


write.csv(data.batch.window.vpa.acute,file="./output_data/vpa_acute_study_window_treatment_data.csv")



