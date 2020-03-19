library(magick)
library(umap)
library(HTSanalyzeR)
library(GSEABase)
library(org.Hs.eg.db)
library(parallel)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(magrittr)
library(reshape2)
library(readxl)
library(scales)
library(gsl)

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


predict.gamma.curve <- function(gene.vec,t.vec,conc.vec,parameter.frame){
    model.pred.points <- as_data_frame(expand.grid(gene=gene.vec,conc=conc.vec,stringsAsFactors=FALSE))
    model.pred.points %>%
        left_join(parameter.frame,by="gene") %>%
        group_by(gene) %>%
        mutate(gc1=create.gamma.curve(conc=conc,gp=g1p[1],g=g1[1],k=k01[1],n=n1[1]),gc2=create.gamma.curve(conc=conc,gp=g2p[1],g=g2[1],k=k02[1],n=n2[1])) %>%
        dplyr::select(gene,conc,gc1,gc2) ->
        model.prediction
}

predict.genes <- function(gene.vec,t.vec,conc.vec,parameter.frame){
    model.pred.points <- as_data_frame(expand.grid(gene=gene.vec,t=t.vec,conc=conc.vec,stringsAsFactors=FALSE))
    model.pred.points %>%
    left_join(parameter.frame,by="gene") %>%
        group_by(gene) %>%
        mutate(model.response=bateman.model(t=t,conc=conc,parms=c(g1p[1],g2p[1],k01[1],n1[1],k02[1],n2[1],a[1],b[1],g1[1],g2[1],n[1],dt[1]))) %>%
        mutate(y=model.response) %>%
        dplyr::select(gene,t,conc,y) ->
        model.prediction
}

collect.kinetic.dose.data <- function(gene.vec,data.fr){
    data.fr=ungroup(data.fr) %>%
        filter(gene %in% gene.vec) %>%
        mutate(symbol=ga[["ensg.symbol"]][gene])

    data.fr %>%
        filter(conc %in% c(0,0.6)) %>%
        group_by(gene,t) %>%
        mutate(y=mny,fc=mny-mny[conc==0]) %>%
        dplyr::select(gene,t,conc,y,fc,sdy) ->
        data.kinetic

    data.fr %>%
        filter(!(conc ==0.6)) %>%
        group_by(gene,t) %>%
        mutate(fc=mny,y=mny+mny[conc==0]) %>%
        dplyr::select(gene,t,conc,y,fc,sdy) %>%
        filter(conc!=0) ->
        data.dose

    data.all <- bind_rows(data.kinetic,data.dose)
}    



combine.data.and.prediction <- function(gene.vec,t.vec,conc.vec,parameter.frame,data.fr){
    data.all <- collect.kinetic.dose.data(gene.vec=gene.vec,data.fr=data.fr)
    model.prediction <- predict.genes(gene.vec,t.vec=t.vec,conc.vec=conc.vec,parameter.frame=parameter.frame)
    model.prediction %>%
        group_by(gene,t) %>%
        mutate(fc=y-y[conc==0]) ->
        model.prediction
        
    model.prediction$sdy <- NA
    model.and.data <- bind_rows(model=model.prediction,data=data.all,.id="type")
}    



mytheme.model.fit.paper <- theme_classic(base_size=10,base_family="Helvetica")+
    theme(legend.key.height = unit(0.75,"lines"))+
    theme(strip.background = element_blank())


load(symbol.conversion.file)

all.genes <- read_csv(all.gene.file,col_names=FALSE) %$%
    X1
regulated.genes <- all.genes[1:n.top.regulated]

colfunc <- colorRampPalette(c("darkblue", "darkorange"))
