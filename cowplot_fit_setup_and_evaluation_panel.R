
args = commandArgs(trailingOnly=TRUE)
data.file=args[1]
all.gene.file=args[2]
symbol.conversion.file=args[3]
opt.file=args[4]
n.top.regulated=as.numeric(args[5])
out.file=args[6]

source("common_setup.R")


######################################################################################################
########### Figure 2 A
###########
######################################################################################################

all.data.batch.corrected <- read_csv(data.file)

gene.symbols.dose <- c("NEFL","SIX3","KLF5","NR2F2","HPGD","GAD1")

ct.genes <- ga[["symbol.ensg"]][gene.symbols.dose]

all.data.batch.corrected %>%
    filter(gene %in% ct.genes) %>%
    filter(batch==3|(conc==0 & t==144)) %>%
    mutate(mny=ifelse(conc==0,0,mny)) %>%
    mutate(symbol=ga[["ensg.symbol"]][gene]) %>%
    mutate(sdy=ifelse(conc==0,0,sdy)) %>%
    mutate(ymin=mny-sdy,ymax=mny+sdy) ->
    plot.object

equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
      yy <- seq(floor(min(x)),ceiling(max(x)),1)
      print(length(yy))
      if (length(yy)<5){
          yy <- seq(floor(min(x)),ceiling(max(x)),1)
      }
      yy
      }
}

plot.object %>%
    mutate(symbol=factor(symbol,levels=gene.symbols.dose,ordered=T)) %>%
    ggplot(aes(x=conc,y=mny,ymin=ymin,ymax=ymax))+
    geom_line(colour="gray40",linetype=1)+
    geom_errorbar(size=0.2,width=0.1)+
    mytheme.model.fit.paper+
    facet_wrap(~symbol,ncol=2,scales="free_y")+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_y_continuous("Fold change at endpoint [Log2]",breaks = pretty_breaks(n = 2))+
    scale_x_continuous("VPA [mM]",breaks=c(0,0.5,1))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(strip.background = element_blank())+
    theme(strip.text = element_text(vjust = -1))  +
    theme(panel.margin.y = unit(0, "lines")) ->
    p.dose




######################################################################################################
########### Figure 2 B
###########
######################################################################################################

gamma.curve.and.model.solution <- function(parameters){
    conc.vec <- c(0,0.2,0.35,0.6)
    par.frame <- as_tibble(t(as.matrix(parameters)))
    par.frame$gene="foo"
    model.prediction <- predict.genes("foo",t.vec=seq(0,144,1),conc.vec=conc.vec,parameter.frame=par.frame)
    gamma.curve.prediction <- predict.gamma.curve("foo",conc.vec=round(seq(0,1,0.01),2),parameter.frame=par.frame)
    list(kinetic.prediction=model.prediction,gamma.curve=gamma.curve.prediction)
}

endpoint.plot <- function(model.prediction.true,model.prediction.wrong){
    min.model <- min(model.prediction.true$y)
    max.model <- max(model.prediction.true$y)+0.2
    margin.vec <- unit(c(0,0,0,0),"pt")
    conc.vec <- unique(model.prediction.true$conc)

    p.schem.ep <- ggplot(model.prediction.wrong,aes(x=t,y=y,colour=as.factor(conc),group=conc))+
        geom_point(data=filter(model.prediction.true,t ==144))+
        geom_line()+
        mytheme.model.fit.paper+
        scale_color_manual("VPA [mM]",values=colfunc(length(conc.vec)))+
        scale_x_continuous("Time")+
        scale_y_continuous("Expression")+
        theme(panel.background = element_blank())+
        theme(legend.position=c(0.2,0.8))+
        theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
        coord_cartesian(ylim=c(min.model,max.model))+
        theme(plot.margin = margin.vec)
}

kinetic.plot <- function(model.prediction.true,model.prediction.wrong){
    margin.vec <- unit(c(0,0,0,0),"pt")
    conc.vec <- unique(model.prediction.true$conc)
    p.schem.kin <- ggplot(filter(model.prediction.wrong,conc %in% c(0,0.35)),aes(x=t,y=y,colour=as.factor(conc),group=conc))+
        geom_point(data=filter(model.prediction.true,conc %in% c(0,0.35),t %in% seq(0,144,24)),aes(shape="data"))+
        geom_line(aes(linetype="model"))+mytheme.model.fit.paper+scale_color_manual("VPA concentration",values=c(colfunc(length(conc.vec))[1],colfunc(length(conc.vec))[which(conc.vec==0.35)]))+
        scale_x_continuous("Time")+
        scale_y_continuous("Expression")+
        scale_linetype_discrete("")+
        scale_shape_discrete("")+
        guides(colour=FALSE)+
        theme(panel.background = element_blank())+
        theme(axis.title.y=element_blank(),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
        theme(plot.margin = margin.vec)+
        theme(legend.background=element_rect(fill="transparent"))
}

gamma.plot <- function(gamma.curve.prediction){
    margin.vec <- unit(c(0,0,0,0),"pt")
    conc.vec <- c(0,0.2,0.35,0.6)
    data.schem=filter(gamma.curve.prediction,conc %in% conc.vec)
    data.schem$ymin=0
    data.schem$ymax=data.schem$gc2
    p.schem.gamma <- ggplot(gamma.curve.prediction,aes(x=conc,y=gc2))+
        geom_line()+
        mytheme.model.fit.paper+
        scale_x_continuous("VPA concentration")+
        labs(y=expression(g[off]~rate))+
        geom_linerange(aes(ymin=ymin,ymax=ymax),data=data.schem)+
        geom_point(colour=colfunc(length(conc.vec)),data=data.schem)+
        theme(panel.background = element_blank())+
        coord_cartesian(xlim=c(-0.15,1))+
        theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
        theme(plot.margin = margin.vec)
}



all.data.batch.corrected <- read_csv(data.file)
opt.par.frame <- read_csv(file=opt.file)

partrue <- c(g1p=6.38,g2p=10,k01=10,n1=7,k02=0.5,n2=3,a=21,b=2e4,g1=0.1,g2=0.29,n=15,dt=0)
parlist <- list(first=partrue,second=partrue,third=partrue)
parlist[["first"]]["n2"] <- 1
parlist[["first"]]["g2p"] <- 6
parlist[["second"]]["k02"] <- 0.55
parlist[["second"]]["n2"] <- 4


prediction.list <- lapply(parlist,gamma.curve.and.model.solution)
kinetic.prediction <- lapply(prediction.list,"[[","kinetic.prediction")
gamma.curve <- lapply(prediction.list,"[[","gamma.curve")

kinetic.plots <- lapply(kinetic.prediction,kinetic.plot,model.prediction.true=kinetic.prediction[["third"]])
dose.plots <- lapply(kinetic.prediction,endpoint.plot,model.prediction.true=kinetic.prediction[["third"]])
gamma.plots <- lapply(gamma.curve,gamma.plot)

all.plots <- lapply(names(parlist),function(x){list(gamma.plots[[x]],dose.plots[[x]],kinetic.plots[[x]])})
all.plots <- do.call("c",all.plots)
all.plots[1:6] <- lapply(all.plots[1:6],function(x){x+theme(axis.title.x=element_blank())})
all.plots[c(2,3,5,6)] <- lapply(all.plots[c(2,3,5,6)],function(x){x+guides(colour=F)})


p.all.cols <- plot_grid(plotlist=all.plots,nrow=3,rel_heights=c(1,1,1.15),rel_widths=c(1,1,1.5))

p.all.cols <- plot_grid(plot_grid(ggdraw()+draw_label(expression(paste("Rate ",'g'['off'],"(VPA)",sep="")),size=10),ggdraw()+draw_label("endpoint data",size=10),ggdraw()+draw_label("kinetic data",size=10),ncol=3),p.all.cols,nrow=2,rel_heights=c(0.1,1))

#p.all.cols <- plot_grid(NULL,p.all.cols,NULL,plot_grid(legend.kin,legend.kin.2,ncol=1),rel_widths=c(0.1,1,0.1,0.2),ncol=4)

######################################################################################################
########### Figure 2 D
###########
######################################################################################################

provide.gamma.fit.data <- function(model.gene){
    model.and.data <- combine.data.and.prediction(gene.vec=ga[["symbol.ensg"]][model.gene],t.vec=seq(0,144,1),conc.vec=round(seq(0,1,0.01),2),parameter.frame=opt.par.frame,data.fr=all.data.batch.corrected)
    gamma.curve.prediction <- predict.gamma.curve(ga[["symbol.ensg"]][model.gene],conc.vec=round(seq(0,1,0.01),2),parameter.frame=opt.par.frame)
    list(model.and.data=model.and.data,gamma.curve=gamma.curve.prediction)
}


kinetic.plot <- function(model.and.data,model.gene){
    ggplot(filter(model.and.data,conc %in% c(0,0.6)),aes(x=t,y=y,ymax=y+sdy,ymin=y-sdy,colour=as.factor(conc),group=conc))+
        geom_point(data=filter(model.and.data,conc %in% c(0,0.6),type=="data"))+
        geom_errorbar(data=filter(model.and.data,conc %in% c(0,0.6),type=="data"))+
        geom_line(data=filter(model.and.data,conc %in% c(0,0.6),type=="model"))+
        annotate(geom="text",x=25,y=Inf,label=model.gene,vjust=1)+
        mytheme.model.fit.paper+
        scale_color_manual("VPA [mM]",values=c(colfunc(4)[1],colfunc(4)[3]))+
        scale_x_continuous("Time [h]")+
        scale_y_continuous("Expression")+
        guides(colour=FALSE)+
        theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) ->
        p.kinetic
}

dose.plot <- function(model.and.data){
    ggplot(filter(model.and.data,t == 144),aes(x=conc,y=fc,ymax=fc+sdy,ymin=fc-sdy))+
        geom_point(data=filter(model.and.data,t == 144,type=="data"))+
        geom_errorbar(data=filter(model.and.data,t == 144,type=="data"))+
        geom_line(data=filter(model.and.data,t == 144,type=="model"))+
        mytheme.model.fit.paper+
        scale_x_continuous("VPA [mM]")+
        scale_y_continuous("Log2 fc")+
        theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) ->
        p.dose
}

g2.plot <- function(gamma.curve.prediction){
    min.curve <- min(gamma.curve.prediction$gc2)
    y.values <- filter(gamma.curve.prediction,conc%in%c(0,0.6))$gc2
    data.schem=tibble(ymax=y.values,ymin=rep(min.curve,2),conc=c(0,0.6),gc2=y.values)
    ggplot(gamma.curve.prediction,aes(x=conc,y=gc2))+
        geom_line()+
        geom_linerange(aes(ymin=ymin,ymax=ymax),data=data.schem)+
        geom_point(colour=c(colfunc(4)[1],colfunc(4)[3]),data=data.schem)+
#        geom_vline(xintercept=c(0,0.6),colour=c(colfunc(4)[1],colfunc(4)[3]),linetype="solid",size=0.8,alpha=1)+
        mytheme.model.fit.paper+
        scale_x_continuous("VPA [mM]")+
        scale_y_continuous(expression('g'['off'])) +
        theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) ->
        p.gamma2
}

g1.plot <- function(gamma.curve.prediction){
    min.curve <- min(gamma.curve.prediction$gc1)
    y.values <- filter(gamma.curve.prediction,conc%in%c(0,0.6))$gc1
    data.schem=tibble(ymax=y.values,ymin=rep(min.curve,2),conc=c(0,0.6),gc1=y.values)
    ggplot(gamma.curve.prediction,aes(x=conc,y=gc1))+
        geom_line()+
        geom_linerange(aes(ymin=ymin,ymax=ymax),data=data.schem)+
        geom_point(colour=c(colfunc(4)[1],colfunc(4)[3]),data=data.schem)+
#        geom_vline(xintercept=c(0,0.6),colour=c(colfunc(4)[1],colfunc(4)[3]),linetype="solid",size=0.8,alpha=1)+
        mytheme.model.fit.paper+
        scale_x_continuous("VPA [mM]")+
        scale_y_continuous(expression('g'['on'])) +
        theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) ->
        p.gamma1
}

gamma.fit.data <- lapply(c("PAX6","LGI1"),provide.gamma.fit.data)
names(gamma.fit.data) <- c("PAX6","LGI1")
model.and.data <- lapply(gamma.fit.data,"[[","model.and.data")
gamma.curve <- lapply(gamma.fit.data,"[[","gamma.curve")

kinetic.plots <- lapply(names(model.and.data),function(x){kinetic.plot(model.and.data[[x]],x)})
names(kinetic.plots) <- names(model.and.data)
dose.plots <- lapply(model.and.data,dose.plot)
g1.plots <- lapply(gamma.curve,g1.plot)
g2.plots <- lapply(gamma.curve,g2.plot)

plot.list <- lapply(c("PAX6","LGI1"),function(x){list(g1.plots[[x]],g2.plots[[x]],kinetic.plots[[x]],dose.plots[[x]])})
plot.list <- do.call("c",plot.list)
plot.list[1:4] <- lapply(plot.list[1:4],function(x){x+scale_x_continuous("")})

p.fit.example <- plot_grid(plotlist=plot.list,nrow=2)



######################################################################################################
########### Figure 2 E
###########
######################################################################################################



parexample <- c(g1p=8.1,g2p=5.9,k01=1.65,n1=4.13,k02=0.49,n2=2.2,a=10.3,b=26e3,g1=0.31,g2=0.97,n=54.3,dt=85.4)
df.0 <- as_data_frame(t(as.matrix(parexample)))
df.0$gene=ga[["symbol.ensg"]]["DAZL"]

model.prediction <- predict.genes(ga[["symbol.ensg"]]["DAZL"],t.vec=seq(0,144,1),conc.vec=round(seq(0,1,0.2),2),parameter.frame=df.0)

model.prediction %>%
    group_by(conc) %>%
    mutate(tmax=t[which.max(y)],ymax=y[which.max(y)])->
    model.prediction

ggplot(filter(model.prediction,conc %in% c(0,0.6)),aes(x=t,y=y,colour=as.factor(conc),group=conc))+
    geom_line(size=1.5)+
    mytheme.model.fit.paper+
    scale_color_manual("VPA [mM]",values=c(colfunc(4)[1],colfunc(4)[3]))+
    geom_hline(linetype="dotted",yintercept=model.prediction$ymax[model.prediction$conc==0][1],colour=colfunc(4)[1])+
    geom_hline(linetype="dotted",yintercept=model.prediction$ymax[model.prediction$conc==0.6][1],colour=colfunc(4)[3])+
    geom_vline(linetype="dotted",xintercept=model.prediction$tmax[model.prediction$conc==0][1],colour=colfunc(4)[1])+
    geom_vline(linetype="dotted",xintercept=model.prediction$tmax[model.prediction$conc==0.6][1],colour=colfunc(4)[3])+
    guides(colour=FALSE)+
    scale_x_continuous("Time")+
    scale_y_continuous("Expression")+
    coord_cartesian(ylim=c(3.5,12))+
    theme(axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank()) ->
    p.vpa.effect.def



######################################################################################################
########### Figure 2 F,G
###########
######################################################################################################


model.prediction <- predict.genes(intersect(regulated.genes,opt.par.frame$gene),t.vec=seq(0,144,1),conc.vec=round(seq(0,1,0.1),2),parameter.frame=opt.par.frame)


model.prediction %>%
    group_by(gene,conc) %>%
    summarise(amplitude=max(y)-min(y),peak=t[which.max(y)]) %>%
    mutate(peak0=peak[conc==0],peak=(peak-peak[conc==0])) %>%
    ungroup %>%
    filter(conc==0.6) %>%
    mutate(peak0=cut(peak0,seq(0,144,24),include.lowest=TRUE)) %>%
    group_by(peak0,conc) %>%
    summarise(md=median(peak),ymax=quantile(peak,0.75),ymin=quantile(peak,0.25)) %>%
    ggplot(aes(x=peak0,group=peak0,y=md,ymax=ymax,ymin=ymin))+
    geom_hline(yintercept=0,colour="red",alpha=0.5,size=1)+
    geom_crossbar(width=0.75,fatten=1.5,alpha=0.5)+
    mytheme.model.fit.paper+
    scale_x_discrete("Peak time untreated during day",labels=as.character(0:5))+
    scale_y_continuous("Time shift [h]")+
    coord_cartesian(ylim=c(-20,40))->
    p.time.shift


model.prediction %>%
    group_by(gene,conc) %>%
    summarise(amplitude=max(y)-min(y),peak=t[which.max(y)]) %>%
    mutate(amplitude=(amplitude-amplitude[conc==0])) %>%
    ungroup %>%
    filter(conc==0.6) %>%
    mutate(amplitude=ifelse(abs(amplitude)>2.5,sign(amplitude)*2.5,amplitude)) %>%
    ggplot(aes(x=amplitude))+geom_histogram(bins=26)+mytheme.model.fit.paper+mytheme.model.fit.paper+scale_x_continuous("Log2 Amplitude shift")+
    coord_cartesian(xlim=c(-2.5,2.5))+
    geom_vline(xintercept=0,colour="red",alpha=0.5,size=1) ->
    p.amplitude.shift

model.prediction %>%
    group_by(gene,conc) %>%
    summarise(amplitude=max(y)-min(y),peak=t[which.max(y)]) %>%
    mutate(amplitude=(amplitude-amplitude[conc==0])) %>%
    ungroup %>%
    filter(conc==0.6) %>%
    ungroup %>%
    summarise(fraction.amplitude.greater.0=sum(amplitude>0)/length(amplitude))


######################################################################################################
########### Collect Panel
###########
######################################################################################################


image_read("./input_data/data_matrix_rotated.svg.png") ->
    data.matrix


p.row.1 <- plot_grid(ggdraw()+draw_image(data.matrix),p.dose,p.all.cols+theme(plot.margin = unit(c(5, 5, 5,10),"pt")),nrow=1,labels=c("A","B","C"),align="v",axis="b",rel_widths=c(0.5,0.5,1))+theme(plot.margin=unit(c(5,5,8,5),"pt"))
p.row.2 <- plot_grid(p.fit.example+theme(plot.margin = unit(c(5, 0, 5,10),"pt")),nrow=1,labels=c("D"))+theme(plot.margin=unit(c(5,5,8,5),"pt"))
p.row.3 <- plot_grid(p.vpa.effect.def+theme(plot.margin = unit(c(5, 30, 5, 20),"pt")),p.amplitude.shift+theme(plot.margin = unit(c(5, 20, 5, 20),"pt")),p.time.shift+theme(plot.margin = unit(c(5,5, 5,20),"pt")),nrow=1,labels=c("E","F","G"),align="h",axis="b")


plot_grid(p.row.1,p.row.2,p.row.3,nrow=3,rel_heights=c(1.25,1,0.8))->
    pp
save_plot(filename=out.file,plot=pp,base_height=5.46/2*2.5,base_aspect_ratio=1.8/2.5*2)

######################################################################################################
########### Fit Statistics
###########
######################################################################################################

fit.data <- as_data_frame(read.csv(data.file,stringsAsFactors=FALSE))


genes.to.test <- intersect(regulated.genes,opt.par.frame$gene)
model.and.data <- combine.data.and.prediction(gene.vec=genes.to.test,t.vec=unique(fit.data$t),conc.vec=unique(fit.data$conc),parameter.frame=opt.par.frame,data.fr=fit.data)

group_by(model.and.data,gene) %>%
    filter(conc ==0) %>%
    summarise(cc=cor(y[type=="model"],y[type=="data"])) ->
    cc.untreated

print("For untreated correlated more than 0.9")
print("out of")
nrow(cc.untreated)
print("percentage")
sum(cc.untreated$cc>0.9)/length(cc.untreated$cc)


model.and.data %>%
    filter(conc ==0.6) %>%    
    filter(t %in% c(24,48,72,96,144)) %>%
    filter(gene %in% gene[type=="data" & abs(fc)>1]) %>%
    group_by(gene) %>%
    summarise(cc=cor(y[type=="model"],y[type=="data"])) %>%
    filter(!is.na(cc)) ->
    cc.treated

print("For VPA 0.6 with |fc|>1 at any time point correlated more than 0.9")
print("out of")
nrow(cc.treated)
print("percentage")
sum(cc.treated$cc>0.9)/length(cc.treated$cc)


model.and.data %>%
    filter(conc!=0.6) %>%    
    filter(t ==144) %>%
    filter(gene %in% gene[type=="data" & abs(fc)>1]) %>%
    group_by(gene) %>%
    summarise(cc=cor(y[type=="model"],y[type=="data"])) %>%
    filter(!is.na(cc)) ->
    cc.conc

print("For endpoint with |fc|>1 correlated more than 0.9")
print("out of")
nrow(cc.conc)
print("percentage")
sum(cc.conc$cc>0.9)/length(cc.conc$cc)


