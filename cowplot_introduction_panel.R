args = commandArgs(trailingOnly=TRUE)
data.file=args[1]
all.gene.file=args[2]
symbol.conversion.file=args[3]
n.top.regulated=as.numeric(args[4])
loess.span=as.numeric(args[5])
out.file=args[6]

source("common_setup.R")

all.data.batch.corrected <- as_data_frame(read.csv(data.file,stringsAsFactors=FALSE))


######################################################################################################
########### Figure 1 B
###########
######################################################################################################


all.data.batch.corrected %>%
    filter(gene %in% regulated.genes) %>%
    filter(conc %in% c(0,0.6)) %>%
    select(gene,t,conc,mny) %>%
    unite("t_conc",t,conc) %>%
    spread(t_conc,mny) %>%
    select(-gene) %>%
    as.matrix ->
    umap.in



custom.config = umap.defaults
custom.config$random_state = 798395
custom.config$n_neighbors= 4
custom.config$n_epochs= 4000

umap.dt=umap(t(umap.in),config=custom.config)
umap.kinetic <- bind_cols(tibble(t_conc=rownames(umap.dt$layout)),as_tibble(umap.dt$layout))

umap.kinetic %>%
    separate(t_conc,into=c("t","conc"),sep="_") %>%
    mutate(t=as.numeric(t)) %>%
    ungroup %>%
    arrange(t) %>%
    ggplot(aes(x=V1,y=V2,group=conc,colour=conc,label=t))+
    geom_path(alpha=0.5)+
    geom_point(size=0.5)+
    geom_text(aes(x=V1+0.1,y=V2+0.1),size=3,show.legend=FALSE)+
    scale_x_continuous("UMAP 1")+
    scale_y_continuous("UMAP 2")+
    scale_colour_manual("VPA [mM]",values=c(colfunc(4)[1],colfunc(4)[3]))+
    mytheme.model.fit.paper+
    theme(legend.position=c(0.3,0.3))->
    p.umap




######################################################################################################
########### Figure 1 C
###########
######################################################################################################

all.data.batch.corrected %>%
    filter(gene %in% regulated.genes) %>%
    filter(conc == 0) %>%
    group_by(gene) %>%
    do(ll=data_frame(t=seq(0,144,0.1),y=predict(loess(mny~t,data=.,span=loess.span),newdata=seq(0,144,0.1)))) %$%
    setNames(ll,gene) ->
    interpolated.data


gene.maximal.reg <- unlist(lapply(interpolated.data,function(x){x$t[which.max(x$y)]}))

all.data.batch.corrected %>%
    filter(gene %in% regulated.genes) %$%
    length(unique(gene)) ->
    no.genes


all.data.batch.corrected %>%
    filter(gene %in% regulated.genes) %>%
    filter(conc %in% c(0,0.6)) %>%
    mutate(gene=factor(gene,levels=names(gene.maximal.reg)[order(gene.maximal.reg,decreasing=TRUE)],ordered=TRUE)) %>%
    group_by(gene) %>%
    mutate(y=(mny-mean(mny))/sd(mny)) ->
    plot.data.hm

plot.data.hm %>%
    filter(conc==0) %>%
    ungroup %>%
    mutate(conc=paste("VPA=",conc,"mM",sep="")) %>%
    ggplot(aes(x=as.factor(t),y=gene,fill=y))+
    geom_tile()+
    scale_fill_gradient2("Z-score\nexpression",low="darkblue",high="darkred")+
    scale_x_discrete("Time [h]")+
    scale_y_discrete(paste("Genes sorted by peak time\n(n=",no.genes,")",sep=""))+
    mytheme.model.fit.paper+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.margin=margin(t=0, r=0, b=0, l=-0.1, unit="cm"))+
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))+
    facet_wrap(~conc)->
    p.heat.1

plot.data.hm %>%
    filter(conc==0.6) %>%
    ungroup %>%
    mutate(conc=paste("VPA=",conc,"mM",sep="")) %>%
    ggplot(aes(x=as.factor(t),y=gene,fill=y))+
    geom_tile()+
    scale_fill_gradient2("Z-score\nexpression",low="darkblue",high="darkred")+
    scale_x_discrete("Time [h]")+
    scale_y_discrete("")+
    mytheme.model.fit.paper+
    theme(axis.title.y=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.margin=margin(t=0, r=0, b=0, l=-0.1, unit="cm"))+
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))+
    facet_wrap(~conc)->
    p.heat.2


## all.data.batch.corrected %>%
## #    filter(gene %in% regulated.genes) %>%
##     filter(gene %in% names(gene.maximal.reg)) %>%
##     filter(t %in% c(24,48,72,96,144)) %>%
##     filter(conc %in% c(0,0.6)) %>%
##     group_by(gene,t) %>%
##     summarise(dy=mny[conc==0.6]-mny[conc==0]) %>%
##     ungroup %>%
##     mutate(gene=factor(gene,levels=names(gene.maximal.reg)[order(gene.maximal.reg,decreasing=TRUE)],ordered=TRUE)) %>% 
##     group_by(gene) %>%
## #    mutate(dy=scale(dy)) %>%
##     mutate(dy=dy/max(abs(dy))) %>%
##     ggplot(aes(x=as.factor(t),y=gene,fill=dy))+
##     geom_tile()+
##     scale_fill_gradient2("Relative\nexpression",low="darkblue",high="darkred")+
##     scale_x_discrete("Time [h]")+
##     scale_y_discrete(paste("Genes (n=",no.genes,")",sep=""))+
##     mytheme.model.fit.paper+
##     theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
##     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
##     theme(legend.margin=margin(t=0, r=0, b=0, l=-0.1, unit="cm"))+
##     guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) ->
##     p1

## save_plot(p1,filename="heatmap_vpa_0.6_fold_change_vs_untreated.pdf")

######################################################################################################
########### Figure 1 D
###########
######################################################################################################

all.data.batch.corrected %>%
    filter(gene %in% regulated.genes) %>%
    filter(t %in% c(24,48,72,96,144)) %>%
    filter(conc %in% c(0,0.6)) %>%
    group_by(gene,t) %>%
    group_by(gene,t) %>%
    summarise(dy=mny[conc==0.6]-mny[conc==0]) %>%
    group_by(gene) %>%
    summarise(dy=max(abs(dy))) %>%
    filter(dy>2) %$%
    gene ->
    strongly.affected.genes

    
gene.maximal.reg <- unlist(lapply(interpolated.data,function(x){x$t[which.max(x$y)]}))
gene.maximal.reg <- gene.maximal.reg[gene.maximal.reg>0 & gene.maximal.reg< 144]


enframe(gene.maximal.reg) %>%
    mutate(strongly.affected=ifelse(name %in% strongly.affected.genes,T,F)) %>%
    mutate(sum.strongly.affected=sum(strongly.affected),sum.not.strongly.affected=sum(!strongly.affected)) %>%
    mutate(strongly.affected.string=ifelse(strongly.affected,paste("Max FC > 2 (n=",sum.strongly.affected[1],")",sep=""),paste("Max FC < 2 (n=",sum.not.strongly.affected[1],")",sep=""))) ->
    plot.object


ggplot(plot.object,aes(x=value,colour=strongly.affected.string,group=strongly.affected.string))+
    stat_ecdf()+
    mytheme.model.fit.paper+
    scale_linetype_discrete("")+
    theme(legend.position="top")+
    scale_x_continuous("Peak Time [h]")+
    scale_y_continuous("Fraction of genes")+
    scale_colour_manual("",values=c("gray50","darkred")) +
    guides(colour=guide_legend(nrow=2,byrow=TRUE))->
    p.fc


######################################################################################################
########### Figure 1 E
###########
######################################################################################################


gene.symbols <- c("DAZL","OLFM3","SHISA2","GPC4","NANOG","BMP5")
gene.symbols <- c("DAZL","OLFM3","SHISA2","HIST1H2AC","NANOG","BMP5")


all.data.batch.corrected %>%
    filter(gene %in% ga[["symbol.ensg"]][gene.symbols])%>%
    mutate(symbol=ga[["ensg.symbol"]][gene]) ->
    annotated.data

filter(annotated.data,t==0) %>%
    mutate(conc=0.6) %>%
    bind_rows(annotated.data) ->
    plot.data

ggplot(filter(plot.data,conc %in% c(0,0.6)),aes(x=t,y=mny,ymin=mny-sdy,ymax=mny+sdy,colour=as.factor(conc),group=conc))+
    geom_line()+
    geom_errorbar()+
    mytheme.model.fit.paper+
    scale_colour_manual("VPA [mM]",values=c(colfunc(4)[1],colfunc(4)[3]))+
    scale_x_continuous("Time [h]")+
    scale_y_continuous("Expression [Log2]")+
    facet_wrap(~symbol,scales="free_y") +
    theme(legend.position="top")+
    theme(legend.margin=margin(t=0, r=0, b=-0.2, l=0, unit="cm")) ->
    p.kinetic





######################################################################################################
########### Figure 1 F
###########
######################################################################################################


solution.equal.rates <- function(t,g1,k){
    y=(g1*t)^k/gamma(k+1)*exp(-g1*t)
    data.frame(t=t,y=y)
}

solution.g1.g2 <- function(t,g1,g2,n){
    if (g1>g2){
        rv <- exp(-g2*t)*t^n*g1^n*hyperg_1F1(a=n,b=1+n,x=-(g1-g2)*t)/gamma(n+1)
    }else{
        rv <- exp(-g1*t)*t^n*g1^n*hyperg_1F1(a=1,b=1+n,x=(g1-g2)*t)/gamma(n+1)
    }
    data.frame(t=t,y=rv)
}


full.multi.chain <- function(k.vec=seq(1,6,1),g1=0.1,g2=0.07){
    k.max <- max(k.vec)+1
# all states but the last decay with rate g1 and have the equal rates solution
    sol.out <- lapply(k.vec,function(x){solution.equal.rates(seq(0,300,0.1),g1,x)})
    names(sol.out) <- k.vec
    plot.object <- melt(sol.out,id.vars=1:2)
    colnames(plot.object) <- c("t","y","k")
# last state has the general solution
    plot.object <- rbind(plot.object,cbind(solution.g1.g2(seq(0,300,0.1),g1,g2,k.max),k=k.max))
    plot.object$type="untreated"
    plot.object <- mutate(plot.object,mark=ifelse(k==k.max,"yes","no"))
    as_data_frame(plot.object)
}

generate.chain.plot <- function(x){
    filter(plot.object,L1==x) %>%
        ggplot(aes(x=t,y=y,group=k,label=k,colour=mark))+
        geom_line()+
        geom_label(data=filter(plot.annot,L1==x),size=3)+
        scale_colour_manual(values=c("darkgray","darkred"))+
        mytheme.model.fit.paper+
        scale_y_continuous("Fraction in State",breaks = pretty_breaks(n = 4))+
        guides(colour = FALSE)+
        coord_cartesian(xlim=c(0,120),ylim=c(0,0.55))+
        theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
        scale_x_continuous("Time")
}


ll <- list()
g1l <- 0.04
g1fact <- 2
g2l <- 0.025
g2fact <- 3

ll[["g1l_g2h"]] <- full.multi.chain(k.vec=0:3,g1=g1l,g2=g2fact*g2l)
ll[["g1h_g2h"]] <- full.multi.chain(k.vec=0:3,g1=g1fact*g1l,g2=g2fact*g2l)
ll[["g1l_g2l"]] <- full.multi.chain(k.vec=0:3,g1=g1l,g2=g2l)
ll[["g1h_g2l"]] <- full.multi.chain(k.vec=0:3,g1=g1fact*g1l,g2=g2l)


plot.object <- bind_rows(ll,.id="L1") %>%
    mutate(L1=factor(L1,levels=c("g1l_g2l","g1h_g2l","g1l_g2h","g1h_g2h"),ordered=TRUE))
                     
margin.vec <- unit(c(1,1,1,1),"mm")

group_by(plot.object,L1,k) %>%
    summarise(t=t[which.max(y)],y=max(y)+0.025,type=type[1],mark=mark[1]) ->
    plot.annot


p.chain <- ggplot(plot.object,aes(x=t,y=y,group=k,label=k,colour=mark))+
    geom_line()+
    geom_label(data=plot.annot,size=3)+
    scale_colour_manual(values=c("darkgray","darkred"))+
    mytheme.model.fit.paper+
    scale_y_continuous("Fraction in State",breaks = pretty_breaks(n = 4))+
    guides(colour = FALSE)+
    coord_cartesian(xlim=c(0,120),ylim=c(0,0.55))+
    facet_wrap(~L1,ncol=2) +
    theme(strip.text = element_blank()) +
    theme(panel.spacing = unit(0.3, "lines"))


idx <- c("g1h_g2l","g1l_g2l","g1l_g2h")
plot.list.chain <- lapply(idx,generate.chain.plot)
plot.list.chain[c(1,3)] <- lapply(plot.list.chain[c(1,3)],function(x){x+theme(axis.title.x=element_blank(),axis.title.y=element_blank())})

image_read("./input_data/g90519.png") -> drawing.chain

p.whole.chain <- plot_grid(ggdraw()+draw_image(drawing.chain,scale=0.98),plot.list.chain[[1]] + annotate("text", x=80, y=0.54, label="Increased~g[on]", parse=TRUE,size=3.6),plot.list.chain[[2]],plot.list.chain[[3]]+annotate("text", x=80, y=0.45, label="Increased~g[off]", parse=TRUE,size=3.6),align="hv")


######################################################################################################
########### Collect Panels
###########
######################################################################################################


image_read("./input_data/Fig1a.jpg") -> drawing.protocol

#p.row.1 <- plot_grid(ggdraw()+draw_image(drawing.protocol),p.umap,p.heat.1+theme(legend.position="none"),p.heat.2,ncol=4,labels=c("A","B","C"),rel_widths=c(1,0.8,0.65,0.5))

p.row.1 <- plot_grid(ggdraw()+draw_image(drawing.protocol),p.umap,p.heat.1+theme(legend.position="none"),p.heat.2,ncol=4,labels=c("A","B","C"),rel_widths=c(1.5,0.8,0.65,0.55))
p.row.2 <- plot_grid(p.fc,p.kinetic,p.whole.chain,align="v",axis="bt",ncol=3,labels=c("D","E","F"),rel_widths=c(0.5,0.8,1))
## plot_grid(p.row.1,p.row.2,align="hv",nrow=2,rel_heights=c(1,1)) ->
##     pp

plot_grid(p.row.1,p.row.2,align="hv",nrow=2,rel_heights=c(0.9,1)) ->
    pp

## p.row.1 <- plot_grid(ggdraw()+draw_image(drawing.protocol),ncol=1,nrow=1,labels=c("A"))
## p.row.2 = plot_grid(p.umap,p.heat.1+theme(legend.position="none"),p.heat.2,ncol=4,labels=c("B","C"),rel_widths=c(0.8,0.65,0.55))
## p.row.3 <- plot_grid(p.fc,p.kinetic,p.whole.chain,align="v",axis="bt",ncol=3,labels=c("D","E","F"),rel_widths=c(0.5,0.8,1))
## plot_grid(p.row.1,p.row.2,p.row.3,align="hv",nrow=3,rel_heights=c(1,1)) ->
##     pp

save_plot(filename=out.file,plot=pp,base_height=5.46,base_aspect_ratio=1.8)
