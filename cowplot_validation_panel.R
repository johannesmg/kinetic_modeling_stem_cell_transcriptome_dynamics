args = commandArgs(trailingOnly=TRUE)
data.file=args[1]
all.gene.file=args[2]
symbol.conversion.file=args[3]
opt.file=args[4]
pcr.data=args[5]
n.top.regulated=as.numeric(args[6])
out.file=args[7]

source("common_setup.R")

######################################################################################################
########### Figure 3 A
###########
######################################################################################################

prepare.data.model.fit <- function(gene.list,optimal.parameter.frame,data.object){
    model.and.data <- combine.data.and.prediction(gene.vec=gene.list,t.vec=seq(0,144,1,),conc.vec=round(seq(0,1,0.01),2),parameter.frame=optimal.parameter.frame,data.fr=data.object)
    
    model.and.data %>%
        group_by(gene) %>%
        mutate(y=y-y[conc==0&t==0 &type=="data"]) %>%
        ungroup ->
        model.and.data

}




validation.kinetic.plot <- function(model.and.data,margin.vec){
    ct.gene=model.and.data$gene[1]
    gene.string=paste(ga[["ensg.symbol"]][ct.gene],"\n","Log2 fc",sep="")
    ggplot(filter(model.and.data,conc %in% c(0,0.6)),aes(x=t,y=y,ymax=y+sdy,ymin=y-sdy,colour=as.factor(conc),group=conc))+
        geom_point(data=filter(model.and.data,conc %in% c(0,0.6),type=="data"))+
        geom_errorbar(data=filter(model.and.data,conc %in% c(0,0.6),type=="data"))+
        geom_line(data=filter(model.and.data,conc %in% c(0,0.6),type=="model"))+
        mytheme.model.fit.paper+
        scale_color_manual("VPA [mM]",values=c(colfunc(4)[1],colfunc(4)[3]))+
        scale_x_continuous("")+
        scale_y_continuous(gene.string)+
        theme(plot.margin = unit(margin.vec, "cm"))+
        labs(x=NULL, y=NULL)+
        theme(panel.background = element_blank()) +
        guides(colour=FALSE)->
        p.kinetic
}

validation.dose.plot <- function(model.and.data,margin.vec){
    ggplot(filter(model.and.data,t == 144),aes(x=conc,y=fc,ymax=fc+sdy,ymin=fc-sdy))+
        geom_point(data=filter(model.and.data,t == 144,conc!=0.6,type=="data"))+
        geom_errorbar(data=filter(model.and.data,t == 144,conc!=0.6,type=="data"))+
        geom_line(data=filter(model.and.data,t == 144,type=="model"))+
        mytheme.model.fit.paper+
        scale_x_continuous("",breaks=c(0,0.5,1))+
        theme(plot.margin = unit(margin.vec, "cm"))+
        labs(x=NULL, y=NULL)+
        theme(panel.background = element_blank()) +
        scale_y_continuous("")->
        p.dose
}
    
validation.prediction.plot <- function(model.prediction,margin.vec){
    model.prediction %>%
        mutate(y=y-y[t==0 & conc==0]) %>%
        ggplot(aes(x=t,y=y,colour=as.factor(conc),group=conc))+
        geom_line()+
        mytheme.model.fit.paper+
        scale_color_manual("VPA [mM]",values=colfunc(4))+
        scale_x_continuous("")+
        scale_y_continuous("")+
        theme(plot.margin = unit(margin.vec, "cm"))+
        labs(x=NULL, y=NULL)+
        theme(panel.background = element_blank()) +
        guides(colour=FALSE)->
        p.pred
}

validation.validation.plot <- function(validation.data,margin.vec){
    ggplot(filter(validation.data,t>=0),aes(x=t,y=y,colour=as.factor(conc),ymax=y+sdy,ymin=y-sdy,group=conc))+
        geom_point(size=0.5)+
        geom_errorbar(size=0.5,width=0.5)+
        geom_line(size=0.5)+
        mytheme.model.fit.paper+
        scale_color_manual("VPA [mM]",values=colfunc(4))+
        scale_x_continuous("")+
        scale_y_continuous("")+
        theme(plot.margin = unit(margin.vec, "cm"))+
        labs(x=NULL, y=NULL)+
        theme(panel.background = element_blank()) +
        guides(colour=FALSE)->
        p.valid
}

arrange.validation.panel <- function(genes,kinetic.plots,dose.plots,prediction.plots,validation.plots){
    panel.plots <- lapply(genes,function(x){list(kinetic.plots[[x]],dose.plots[[x]],prediction.plots[[x]],validation.plots[[x]])})
    panel.plots <- do.call("c",panel.plots)
    panel.plots[1:20] <- lapply(panel.plots[1:20],function(x){x+theme(axis.text.x=element_blank(),axis.title.x=element_blank())})
    panel.plots[c(21,23,24)] <- lapply(panel.plots[c(21,23,24)],function(x){x+scale_x_continuous("Time [h]")})
    panel.plots[[22]] <- panel.plots[[22]]+scale_x_continuous("VPA [mM]",breaks=c(0,0.5,1))
    panel.plots <- lapply(panel.plots,function(x){x+theme(plot.background=element_blank(),panel.background=element_blank())})
    plot_grid(plotlist=panel.plots,ncol=4,align="vh",axis="btlr")
}


opt.par.frame <- read_csv(file=opt.file)
fit.data <- read_csv(data.file)
all.data.batch.corrected <- fit.data


pcr.data.utr <- read_excel(path=pcr.data,sheet=1)
pcr.data.utr <- gather(pcr.data.utr,key="replicate",value="dct",4:6)
translation.vector <- c(0.35,0.6,0.8,0)
names(translation.vector) <- names(table(pcr.data.utr$condition))
pcr.data.utr$condition <- translation.vector[pcr.data.utr$condition]
pcr.data.utr <- na.omit(pcr.data.utr)
pcr.data.utr$DoD <- as.numeric(sapply(pcr.data.utr$DoD,substr,start=2,stop=5))

group_by(pcr.data.utr,Gene,replicate) %>%
    mutate(dct=-dct+dct[DoD==-3]) ->
    pcr.data.utr

mean.pcr.data <- group_by(pcr.data.utr,Gene,DoD,condition) %>%
    summarise(sddct=sd(dct,na.rm=TRUE),dct=mean(dct,na.rm=TRUE))

mean.pcr.data$t <- mean.pcr.data$DoD*24
mean.pcr.data$y <- mean.pcr.data$dct

day.zero <- as.data.frame(expand.grid(Gene=unique(mean.pcr.data$Gene),DoD=0,condition=unique(mean.pcr.data$condition),dct=0,t=0,y=0,sddct=0))
mean.pcr.data <- bind_rows(mean.pcr.data,day.zero)

mean.pcr.data <- mean.pcr.data[,c("t","y","sddct","Gene","condition")]
colnames(mean.pcr.data) <- c("t","y","sdy","gene","conc")
mean.pcr.data$symbol=mean.pcr.data$gene
mean.pcr.data$gene=ga[["symbol.ensg"]][mean.pcr.data$gene]

mean.pcr.data %>%
    filter(!is.na(gene)) ->
    mean.pcr.data

gene.symbols <- c("PAX6","DAZL","OTX2","BMP5","SOX9","POU4F1","GAD1","NPY","OLFM3","PAX3","PIPOX","SNAI2")
gene.symbols <- gene.symbols[ga[["symbol.ensg"]][gene.symbols] %in% opt.par.frame$gene]


model.prediction <- predict.genes(gene.vec=ga[["symbol.ensg"]][gene.symbols],t.vec=seq(0,144,1),conc.vec=round(c(0,0.35,0.6,0.8),2),parameter.frame=opt.par.frame) %>%
    ungroup
data.model.fit <- prepare.data.model.fit(ga[["symbol.ensg"]][gene.symbols],optimal.parameter.frame=opt.par.frame,data.object=all.data.batch.corrected)

group_split(data.model.fit,gene) ->
    data.model.fit.list

names(data.model.fit.list) <- group_keys(data.model.fit,gene)$gene

group_split(model.prediction,gene) ->
    model.prediction.list

names(model.prediction.list) <- group_keys(model.prediction,gene)$gene

group_split(mean.pcr.data,gene) ->
    validation.data.list

names(validation.data.list) <- group_keys(mean.pcr.data,gene)$gene


margin.vector <- rep(-0.1,4)
margin.vector[2] <- -0.6
margin.vector[1] <- -0.2

kinetic.plots <- lapply(data.model.fit.list,validation.kinetic.plot,margin.vec=margin.vector)
dose.plots <- lapply(data.model.fit.list,validation.dose.plot,margin.vec=margin.vector)
prediction.plots <- lapply(model.prediction.list,validation.prediction.plot,margin.vec=margin.vector)
validation.plots <- lapply(validation.data.list,validation.validation.plot,margin.vec=margin.vector)

first.panel.genes <- ga[["symbol.ensg"]][gene.symbols][1:6]
second.panel.genes <- ga[["symbol.ensg"]][gene.symbols][7:12]


row1.col1 <- arrange.validation.panel(first.panel.genes,kinetic.plots,dose.plots,prediction.plots,validation.plots)
row1.col2 <- arrange.validation.panel(second.panel.genes,kinetic.plots,dose.plots,prediction.plots,validation.plots)


mean.pcr.data %>%
    filter(symbol %in% gene.symbols) %>%
    group_by(gene) %>%
    mutate(y=y-y[t==0 & conc==0]) %>%
    rename(data=y) %>%
    left_join(group_by(model.prediction,gene) %>% mutate(prediction=y-y[t==0 & conc==0]),by=c("gene","t","conc")) %>%
    filter(t!=-72) %>%
    group_by(gene,symbol) %>%
    summarise(cr=cor(data,prediction)) %>%
    filter(!is.na(cr)) %>%
    ungroup %>%
    summarise(mean(cr))

ggplot(tibble(t=1:4,y=1:4,cl=c(0,0.35,0.6,0.8)),aes(x=t,y=y,colour=as.factor(cl)))+
    geom_line()+
    geom_point()+
    scale_colour_manual("VPA [mM]", values=colfunc(4))+mytheme.model.fit.paper ->
    p1

save_plot(p1,filename=paste(out.file,"_legend.pdf",sep=""))


######################################################################################################
########### Figure 3 B, C
###########
######################################################################################################



ct.gene <- ga[["symbol.ensg"]]["GAD1"]
model.prediction <- predict.genes(gene.vec=ct.gene,t.vec=seq(0,300,1),conc.vec=round(seq(0,1,0.01),2),parameter.frame=opt.par.frame)
gamma.curve.prediction <- predict.gamma.curve(gene.vec=ct.gene,conc.vec=round(seq(0,1,0.01),2),parameter.frame=opt.par.frame)

no.colors <- length(unique(model.prediction$conc))

model.prediction %>%
    filter(conc %in% seq(0,1,0.08)) %>%
    ggplot(aes(x=t,y=y,colour=conc,group=conc))+
    geom_line()+
    mytheme.model.fit.paper+
    scale_color_gradient("VPA [mM]",low="darkblue",high="darkorange")+
    scale_y_continuous("Expression")+
    scale_x_continuous("Time [h]")+
    geom_vline(xintercept=c(72,96,144),linetype=3)+
    annotate("text",x=144+11,y=10.3,label="d6")+
    annotate("text",x=96+11,y=10.3,label="d4")+
    annotate("text",x=72+11,y=10.3,label="d3")+
    guides(colour = guide_colourbar(barwidth = 0.5))+
    theme(panel.background = element_blank()) +
    coord_cartesian(ylim=c(5,10.8))->
    p.kinetic


model.prediction %>%
    filter(t%in%c(72,96,144)) %>%
    mutate(fc=y-y[conc==0]) %>%
    ggplot(aes(x=conc,y=fc,group=t,linetype=as.factor(t)))+
    geom_path()+
    annotate("text",x=0.8,y=3.4,label="d3")+
    annotate("text",x=0.8,y=2.65,label="d4")+
    annotate("text",x=0.8,y=1,label="d6")+
    mytheme.model.fit.paper+
    scale_color_gradient("VPA [mM]",low="darkblue",high="darkorange")+
    scale_x_continuous("VPA [mM]")+
    scale_linetype_manual("Measurement",labels=c("d3","d4","d6"),values=rep(1,3))+
    guides(colour = guide_colourbar(barwidth = 0.5),linetype=FALSE)+
    scale_y_continuous("Log2 fc at day 6")->
    p.fold.change


equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(round(min(x)+d,2), round(max(x)-d,2), length=n)
  }
}


gamma.curve.prediction %>%
    gather("curve","rate",gc1:gc2) %>%
    mutate(label="") %>%
    mutate(label=ifelse(curve=="gc1" & conc==0.5,"gon",label)) %>%
    mutate(label=ifelse(curve=="gc2" & conc==0.5,"goff",label)) %>%
    group_by(curve) %>%
    mutate(dr=max(rate)-(max(rate)-min(rate))/10)->
    rate.plot.in


ggplot(rate.plot.in,aes(x=conc,y=rate,label=label))+
    geom_line()+
    mytheme.model.fit.paper+
    scale_x_continuous("VPA [mM]")+
    scale_y_continuous("Rate [1/s]",breaks = pretty_breaks(n = 3))+
    facet_wrap(~curve,ncol=1,scales="free_y")+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),strip.background = element_blank(), strip.text.x = element_blank())->
    p.rate


lab1 <- expression('g'['on'])
ann_text <- filter(rate.plot.in,conc==0.5,curve=="gc1")

p.rate <- p.rate+
    geom_text(data = ann_text,aes(y=dr),label = lab1)

lab2 <- expression('g'['off'])
ann_text <- filter(rate.plot.in,conc==0.5,curve=="gc2")

p.rate <- p.rate+
    geom_text(data = ann_text,aes(y=dr),label = lab2)


######################################################################################################
########### Collect Panel
###########
######################################################################################################

row.1 <- plot_grid(row1.col1+theme(plot.margin = unit(c(2, 2, 0, 0.1), "cm")),row1.col2+theme(plot.margin = unit(c(2, 2, 0, 0.1), "cm")),ncol=2,align="hv",labels=c("A",""))

row.2 <- plot_grid(p.rate+theme(plot.margin=unit(c(5,20,5,15),"pt")),p.kinetic,p.fold.change,nrow=1,rel_widths=c(0.85,1.5,0.85),labels=c("B","C",""),scale=0.9)
plot_grid(row.1,row.2,nrow=2,rel_heights=c(1,0.25)) ->
    pp

#print(out.file)
save_plot(filename=out.file,plot=pp,base_height=5.46/2*3,base_aspect_ratio=1.8/3*2)
