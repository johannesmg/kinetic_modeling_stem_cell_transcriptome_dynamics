bateman.model <- function(t,conc,parms){
    g1p <- parms[1]
    g2p <- parms[2]
    g1 <- parms[9]
    g2 <- parms[10]
    k01 <- parms[3]
    n1 <- parms[4]
    k02 <- parms[5]
    n2 <- parms[6]
    gg1 <- create.gamma.curve(conc,g1p,g1,k01,n1)
    gg2 <- create.gamma.curve(conc,g2p,g2,k02,n2)
    model.extendedgamma.new(t,gg1,gg2,parms[7:12])
}


model.extendedgamma.new <- function(t,g1,g2,parms){
    A=parms[1]
    B=parms[2]
    n=parms[5]
    t=t+parms[6]
    #make use of Kummer's transformation to prevent
                                        #machine size infinities
    g1.greater <- g1>g2
    g1.smaller <- !g1.greater
    rv <- rep(0,length(g1))
    if (sum(g1.greater)>0){
        rv[g1.greater] <- A+B*exp(-g2[g1.greater]*t[g1.greater])*t[g1.greater]^n*g1[g1.greater]^n*hyperg_1F1(a=n,b=1+n,x=-(g1[g1.greater]-g2[g1.greater])*t[g1.greater])/gamma(n+1)
    }
    if (sum(g1.smaller)>0){
        rv[g1.smaller] <- A+B*exp(-g1[g1.smaller]*t[g1.smaller])*t[g1.smaller]^n*g1[g1.smaller]^n*hyperg_1F1(a=1,b=1+n,x=(g1[g1.smaller]-g2[g1.smaller])*t[g1.smaller])/gamma(n+1)
    }
    log2(rv+1)
}



create.gamma.curve <- function(conc,gp,g,k,n){
                                        # this should only happen for large n
                                        # and large n's should not be taken...
    if (k^n>0 & max(conc)>0){
        rv <- (gp-g)*conc^n/(k^n+conc^n)+g
    }else{
        rv <- rep(g,length(conc))
        names(rv) <- names(conc)
    }
    rv
}

library(gsl)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)


args = commandArgs(trailingOnly=TRUE)
opt.file=args[1]
outfile=args[2]
data.file=args[3]
gene.list.file=args[4]
no.top.regulated=args[5]

read.csv(file=data.file,stringsAsFactors=FALSE) %>%
    as_data_frame ->
    all.data

gene.out.frame=read_csv(file=opt.file)
gene.list=read_csv(file=gene.list.file,col_names=FALSE) %$% X1

gene.list=gene.list[1:no.top.regulated]


present.genes <- intersect(all.data$gene,unique(gene.out.frame$gene))
present.genes <- intersect(present.genes,gene.list)

print("Following genes not found")
print(setdiff(gene.list,present.genes))

model.pred.points <- expand.grid(gene=present.genes,t=seq(0,144,4),conc=seq(0,0.6,0.025),stringsAsFactors=FALSE)

model.pred.points %>%
    left_join(gene.out.frame,by="gene") %>%
    group_by(gene) %>%
    mutate(model.response=bateman.model(t=t,conc=conc,parms=c(g1p[1],g2p[1],k01[1],n1[1],k02[1],n2[1],a[1],b[1],g1[1],g2[1],n[1],dt[1]))) %>%
    group_by(gene,t) %>%
    mutate(fc=model.response-model.response[conc==0]) %>%
    filter(conc!=0) %>%
    select(gene,t,conc,fc) ->
    fc.prediction

all.data %>%
    filter(gene %in% present.genes) %>%
    filter(conc %in% c(0,0.6)) %>%
    group_by(gene,t) %>%
    mutate(fc=mny-mny[conc==0]) %>%
    filter(conc!=0) %>%
    select(gene,t,conc,fc) ->
    fc.data.kinetic

all.data %>%
    filter(gene %in% present.genes) %>%
    filter(!(conc %in% c(0,0.6))) %>%
    group_by(gene,t) %>%
    mutate(fc=mny) %>%
    filter(conc!=0) %>%
    select(gene,t,conc,fc) ->
    fc.data.dose

fc.data <- bind_rows(fc.data.kinetic,fc.data.dose)



bind_rows(data=fc.data,model=fc.prediction,.id="type") %>%
    unite("type_t_conc",type,t,conc) ->
    gsea.data.in

colnames(gsea.data.in) <- c("sample","gene","value")

#write_csv(gsea.data.in,path="./out_lhs_07/gsea_in_model_data.csv")
write_csv(gsea.data.in,path=outfile)
