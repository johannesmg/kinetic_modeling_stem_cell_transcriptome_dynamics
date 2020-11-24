
args = commandArgs(trailingOnly=TRUE)
data.file=args[1]
all.gene.file=args[2]
symbol.conversion.file=args[3]
opt.file=args[4]
full.msigdb=args[5]
gsea.full.msigdb=args[6]
gsea.developmental=args[7]
vpa.window.data=args[8]
vpa.acute.study.data=args[9]
mica.gene.sets=args[10]
kegg.msigdb=args[11]
gsea.out.nc=args[12]
n.top.regulated=as.numeric(args[13])
out.file=args[14]
sig.thresh.gsea.kegg=as.numeric(args[15])
neural.crest.gmt.file=args[16]

source("common_setup.R")

######################################################################################################
########### Figure 4 A
###########
######################################################################################################

replace.terms <- function(term){
    term %>%
        gsub("_"," ",x=.) %>%
        gsub("NEURAL CREST STEM CELL UP","Neural Crest Genes\n",x=.) %>%
        gsub("LEE","",x=.) ->
        rv
    rv
     }


gsea.all <- read_csv(file=gsea.full.msigdb,col_names=F)
colnames(gsea.all) <- c("fdr","pvals","score","term","sample","min.set.size","no.reps","set.definition.file")

gsea.all %>%
    separate(sample,into=c("type","t","conc"),sep="_") %>%
    mutate(t=as.numeric(t)) %>%
    mutate(conc=as.numeric(conc)) %>%
    filter(t>0) %>%
    filter(conc>0) %>%
    filter(conc<0.65) %>%
    filter(t!=144|type=="data") %>%
    filter(conc!=0.6|type=="data") %>%
    filter(term%in%term[fdr==min(fdr)]) %>%
    mutate(term=replace.terms(term)) %>%
    complete(t,conc,term) ->
gsea.plot.object



p.gsea.all <- ggplot(gsea.plot.object,aes(x=t,y=conc,fill=1-fdr))+
    geom_tile()+
    scale_fill_gradient2("Significance",low="darkblue",mid="darkblue",midpoint=0.5,high="gold",guide="colourbar",na.value="gray80",breaks=c(min(1-gsea.all$fdr),max((1-gsea.all$fdr))),labels=c("low","high"))+
    scale_x_continuous("Time [h]")+
    scale_y_continuous("VPA [mM]")+
    facet_wrap(~term,nrow=1)+
    mytheme.model.fit.paper+
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5,title.position="left"))+
    theme(legend.title=element_text(angle=90,vjust=0.5,hjust=0.5),legend.title.align=+0.5)+
    theme(legend.margin=margin(l = 0, unit='pt'))


######################################################################################################
########### Figure 4 B
###########
######################################################################################################

replace.terms <- function(term){
    term %>%
        gsub("_"," ",x=.) %>%
        gsub("PATHWAY","",x=.) %>%
        gsub("KEGG ","",x=.) %>%
        gsub("BETA ","BETA\n",x=.) %>%
        gsub("WNT ","WNT\n",x=.) %>%
        gsub("LIGAND ","LIGAND\n",x=.) %>%
        gsub("MAPK ","MAPK\n",x=.) ->
        rv
    rv
     }

gsea.sigdev <- read_csv(file=gsea.developmental,col_names=F)
colnames(gsea.sigdev) <- c("fdr","pvals","score","term","sample","min.set.size","no.reps","set.definition.file")

print("Number of terms")
print(gsea.sigdev %>% distinct(term) %>% nrow)

mutate(gsea.sigdev,dir.enr=sign(score)*(1-pvals)) %>%
    separate(sample,into=c("type","t","conc"),sep="_") %>%
    mutate(t=as.numeric(t)) %>%
    mutate(conc=as.numeric(conc)) %>%
    filter(t>0) %>%
    mutate(term=replace.terms(term)) %>%
    filter(conc>0) %>%
    filter(conc<0.65) %>%
    filter(t!=144|type=="data") %>%
    filter(conc!=0.6|type=="data") %>%
    filter(term%in%term[fdr<sig.thresh.gsea.kegg]) %>%
    complete(t,conc,term) ->
    gsea.plot.object

print(head(gsea.plot.object))

p.gsea.dev <- ggplot(gsea.plot.object,aes(x=t,y=conc,fill=score))+
    geom_tile()+
    scale_fill_gradient2("Enrichment score",low="darkblue",high="darkred",midpoint=0,mid="white",guide="colourbar",na.value="gray80")+
    scale_x_continuous("Time [h]")+
    scale_y_continuous("VPA [mM]")+
    facet_wrap(~term,nrow=2)+
    theme(panel.spacing = unit(1, "lines"))+
    mytheme.model.fit.paper+
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5,title.position="left"))+
    theme(legend.title=element_text(angle=90,vjust=0.5,hjust=0.5),legend.title.align=+1)+
    theme(legend.margin=margin(l = 0, unit='pt'))
   

######################################################################################################
########### Figure 4 D,E
###########
######################################################################################################


perform.wilcox <- function(data.in){
    wilcox.test(x=data.in$fc,alternative="greater")$p.value
}

p.val.to.star <- function(value){
    ifelse(value<0.05,ifelse(value<0.01,ifelse(value<0.001,"***","**"),"*"),"")
}



x <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(x)
ensembl.entrez.mapping <- as.list(x[mapped_genes])
ensembl.entrez.mapping <- melt(ensembl.entrez.mapping)
colnames(ensembl.entrez.mapping) <- c("ensembl","entrez")
ensembl.entrez.mapping[,1] <- as.character(ensembl.entrez.mapping[,1])
ensembl.entrez.mapping[,2] <- as.character(ensembl.entrez.mapping[,2])


gmt.object <- getGmt(full.msigdb)
msigdb.list <- geneIds(gmt.object)

gmt.object <- getGmt(kegg.msigdb)
kegg.list <- geneIds(gmt.object)


all.data.batch.corrected <-  read_csv(data.file)

load(file=mica.gene.sets)

genes.wnt.up=gene.sets[["Wnt_up"]]
genes.wnt.down=gene.sets[["Wnt_dn"]]
neural.crest.ensg <- unique(ensembl.entrez.mapping[ensembl.entrez.mapping$entrez %in% msigdb.list[["LEE_NEURAL_CREST_STEM_CELL_UP"]],"ensembl"])
neural.crest.ensg.dn <- unique(ensembl.entrez.mapping[ensembl.entrez.mapping$entrez %in% msigdb.list[["LEE_NEURAL_CREST_STEM_CELL_DN"]],"ensembl"])
kegg.tgf_beta<- unique(ensembl.entrez.mapping[ensembl.entrez.mapping$entrez %in% kegg.list[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]],"ensembl"])

gene.set.wnt <- list("up/induced"=genes.wnt.up,"down/repressed"=genes.wnt.down,other=setdiff(regulated.genes,union(genes.wnt.up,genes.wnt.down)))
gene.set.tgfb <- list("up/induced"=kegg.tgf_beta,"down/repressed"=c(),other=setdiff(regulated.genes,kegg.tgf_beta))
gene.set.nc <- list("up/induced"=neural.crest.ensg,"down/repressed"=neural.crest.ensg.dn,other=setdiff(regulated.genes,union(neural.crest.ensg,neural.crest.ensg.dn)))

gene.sets <- list("Wnt"=gene.set.wnt,"TGFb"=gene.set.tgfb[c("up/induced","other")],"Neural Crest"=gene.set.nc)
mm <- melt(gene.sets)
colnames(mm) <- c("gene","direction","set")

kinetic.dose.data <- collect.kinetic.dose.data(regulated.genes,data.fr=all.data.batch.corrected)
kinetic.dose.data <- left_join(kinetic.dose.data,mm,by="gene")
kinetic.dose.data <- mutate(kinetic.dose.data,direction=factor(direction,levels=c("up/induced","down/repressed","other"),ordered=TRUE))

plot.annotation.data <- filter(kinetic.dose.data,t%in%c(96,144),conc==0.6,direction=="up/induced") %>%
    group_by(set,direction,conc) %>%
    summarise(mfc=mean(fc))


kinetic.dose.data %>%
    filter(conc==0.6|t==0) %>%
    filter(direction!="other") %>% filter(direction!="down/repressed") %>%
    group_by(t,set,direction) %>%
    do(bind_cols(as_tibble(perform.wilcox(.)),mfc=mean(.$fc))) %>%
    mutate(is.sig=value<0.05) %>%
    mutate(asterisk=p.val.to.star(value)) ->
    is.sig.kinetic

group_by(kinetic.dose.data,set,direction,conc,t) %>%
    summarise(mfc=mean(fc)) %>%
    filter(conc==0.6|t==0) %>%
    filter(direction!="other") %>%
    filter(direction!="down/repressed") %>%
    ggplot(aes(x=t,y=mfc,colour=set,linetype=direction,group=interaction(set,direction)))+
    geom_hline(yintercept=0,colour="gray80")+
    geom_point()+
    geom_line()+
    geom_text(data=plot.annotation.data,size=2.85,aes(x=112,y=mfc+0.22,label=set))+
    geom_text(data=is.sig.kinetic,aes(label=asterisk,y=mfc+ifelse(set=="Wnt",-2.0,ifelse(set=="TGFb",-1.25,1))*(0.1+as.numeric(t<50)*0.05)),angle=90,vjust=0.75,size=4,hjust=0)+
    scale_colour_brewer(palette="Set1")+
    scale_linetype_manual("",values=c("up/induced"="solid","down/repressed"="dashed")) +
    scale_x_continuous("Time [h]")+
    scale_y_continuous("Fc VPA to untreated [Log2]")+
    guides(colour=FALSE)+
    guides(linetype=FALSE)+
    mytheme.model.fit.paper+
    theme(legend.background=element_blank())+
    theme(legend.position=c(0.3,0.16))+
    coord_cartesian(ylim=c(-0.6,1.65))->
    p.kinetic.pathways


plot.annotation.data <- filter(kinetic.dose.data,t==144,conc==0.8,direction=="up/induced") %>%
    group_by(set,direction,conc) %>%
    summarise(mfc=mean(fc))

plot.annotation.data[plot.annotation.data$set=="Neural Crest",]$mfc <- plot.annotation.data[plot.annotation.data$set=="Neural Crest",]$mfc-0.6


kinetic.dose.data %>%
    filter(conc>0) %>%
    filter(conc!=0.6,t==144) %>%
    filter(direction=="up/induced") %>%
    group_by(conc,set,direction) %>%
    do(bind_cols(as_tibble(perform.wilcox(.)),mfc=mean(.$fc))) %>%
    mutate(is.sig=value<0.05) %>%
    mutate(asterisk=p.val.to.star(value)) ->
    is.sig.conc


group_by(kinetic.dose.data,set,direction,conc,t) %>%
    summarise(mfc=mean(fc)) %>%
    filter(t==144,conc!=0.6) %>%
    filter(direction!="other") %>%
    filter(direction!="down/repressed") %>%
    ggplot(aes(x=conc,y=mfc,colour=set,linetype=direction,group=interaction(set,direction)))+
    geom_hline(yintercept=0,colour="gray80")+
    geom_point()+
    geom_line()+
    geom_text(data=plot.annotation.data,size=2.85,aes(x=conc,y=mfc+0.2,label=set))+
    geom_text(data=is.sig.conc,aes(label=asterisk,y=mfc+ifelse(set=="Wnt",-7,1)*(0.04+(1-conc)*0)),angle=90,vjust=0.75,size=4,hjust=0)+
    scale_colour_brewer(palette="Set1")+
    scale_linetype_manual("",values=c("up/induced"="solid","down/repressed"="dashed")) +
    scale_x_continuous("VPA [mM]")+
    scale_y_continuous("Fc VPA to untreated [Log2]")+
    mytheme.model.fit.paper+
    guides(colour=FALSE,linetype=FALSE)+
    coord_cartesian(ylim=c(-0.55,1.65))->
    p.dose.pathways


######################################################################################################
########### Figure 4 F
###########
######################################################################################################


all.data.batch.corrected %>%
    filter(gene %in% regulated.genes) %>%
    filter(t %in% c(48,72,96)) %>%
    filter(conc %in% c(0,0.6)) %>%
    group_by(gene,t) %>%
    summarise(dy=mny[conc==0.6]-mny[conc==0]) %>%
    group_by(gene) ->
    fc.data.kinetic

load(file=mica.gene.sets)

colnames(gene.sets[["Wnt_frame"]]) <- c("gene","z_value","symbol")

fc.data.kinetic %>%
    left_join(gene.sets[["Wnt_frame"]]) %>%
    filter(!is.na(z_value)) %>%
    spread(t,dy) %>%
    mutate(sgn=factor(c("Wnt\nrepressed","Wnt\ninduced")[as.numeric(z_value>0)+1],levels=c("Wnt\nrepressed","Wnt\ninduced"),ordered=TRUE)) %>%
    group_by(sgn) %>%
    ungroup %>%
    filter(abs(z_value)>3.5) %>%
    arrange(z_value) %>%
    gather("t","dy",'48':'96') %>%
    mutate(symbol=factor(symbol,levels=unique(symbol),ordered=T)) %>%
    mutate(t=factor(t,levels=c("48","72","96","120"),ordered=T))->
    plot.data.fc.kinetic

plot.data.fc.kinetic <- bind_rows(plot.data.fc.kinetic,tibble(gene="ENSG00000138083",z_value=NA,symbol=factor("SIX3",levels=unique(plot.data.fc.kinetic$symbol),ordered=T),sgn="Wnt\nrepressed",t=factor("120",levels=c("48","72","96","120"),ordered=T),dy=NA))
divider.line.y.lower <- sum(unique(plot.data.fc.kinetic$z_value)<0,na.rm=T)
divider.line.y.upper <- sum(unique(plot.data.fc.kinetic$z_value)>0,na.rm=T)


plot.data.fc.kinetic %>%
    mutate(sgn.ref=z_value>0,sgn.data=dy>0) %>%
    filter(!is.na(sgn.ref) & !(is.na(sgn.data))) %>%
    dplyr::select(t,sgn.ref,sgn.data) %>%
    group_by(t) %>%
    do(as_tibble(fisher.test(x=xtabs(data=.[,c("sgn.ref","sgn.data")]),alternative="greater")$p.value)) %>%
    mutate(label.sig=paste(t,"\n",p.val.to.star(value),sep="")) ->
    p.val.annotation.wnt.heatmap


plot.data.fc.kinetic %>%
    mutate(sgn.ref=z_value>0,sgn.data=dy>0) %>%
    filter(!is.na(sgn.ref) & !(is.na(sgn.data))) %>%
    dplyr::select(t,sgn.ref,sgn.data) %>%
    group_by(t) %>%
    do(as_tibble(xtabs(data=.[,c("sgn.ref","sgn.data")])))

    

ggplot(plot.data.fc.kinetic,aes(x=t,y=symbol,fill=dy))+
    geom_tile()+
    geom_hline(yintercept=divider.line.y.lower+0.5)+
    scale_fill_gradient2("Log2 fc",low="darkblue",high="darkred",na.value="white")+
    mytheme.model.fit.paper+
    scale_y_discrete("")+
    scale_x_discrete("Time [h]",breaks=c("48","72","96"),labels=p.val.annotation.wnt.heatmap$label.sig)+
    theme(plot.margin = unit(c(5, 0, 5, 5), "pt"))+
    theme(axis.text.y=element_text(size=7))+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))->
    p1


ggplot(plot.data.fc.kinetic,aes(x=as.factor(0),y=factor(sgn),fill=1))+
    geom_tile()+
    scale_fill_gradient2(low="darkblue",high="darkred")+
    mytheme.model.fit.paper+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.line.y=element_blank())+
    scale_y_discrete("",position="right")+
    guides(fill=FALSE) +
    theme(plot.background=element_blank())+
    theme(axis.text.y.right=element_text(angle=270,hjust=0.5,vjust=0,size=10,lineheight=0.75))->
    p2

legend <- get_legend(p1)
prow <- plot_grid(p1+theme(legend.position="none"),p2,align="h",rel_widths=c(1,0.3))
p.wnt.heatmap <- plot_grid( prow, legend, rel_widths = c(1, .25),ncol=2)
p.wnt.heatmap <- p1+ annotate("text", x = 4, y = divider.line.y.lower/2, label = "Wnt repressed",angle=90,hjust=0.5)+ annotate("text", x = 4, y = divider.line.y.lower+divider.line.y.upper/2+0.5, label = "Wnt induced",hjust=0.5,angle=90)


######################################################################################################
########### Supplemental Figure TGFb
###########
######################################################################################################


all.data.batch.corrected %>%
    filter(gene %in% regulated.genes) %>%
    filter(t %in% t[conc==0.6]) %>%
    filter(conc %in% c(0,0.6)) %>%
    group_by(gene,t) %>%
    summarise(dy=mny[conc==0.6]-mny[conc==0]) %>%
    group_by(gene) ->
    fc.data.kinetic


fc.data.kinetic %>%
    mutate(symbol=ga[["ensg.symbol"]][as.character(gene)]) %>%
    filter(gene %in% gene.set.tgfb[["up/induced"]])%>%
    mutate(t=as.factor(t)) %>%
    ggplot(aes(x=t,y=symbol,fill=dy))+
    geom_tile()+
    scale_fill_gradient2("Log2 fc",low="darkblue",high="darkred",na.value="white")+
    mytheme.model.fit.paper+
    scale_y_discrete("")+
    scale_x_discrete("Time [h]") ->
    p1

save_plot(p1,filename="./figures/supplemental_figure_tgfb.pdf")

######################################################################################################
########### Figure 4 D
###########
######################################################################################################


provide.gsea.plot <- function(current.index){
    hits <- names(gsea.list.in[[current.index]])[gsea.list.in[[current.index]] >quantile(gsea.list.in[[current.index]],0.5)]
    gsca <- new("GSCA", listOfGeneSetCollections=list(set=gene.set.list),geneList=gsea.list.in[[current.index]], hits=hits)
    gsca <- preprocess(gsca, species="Hs", initialIDs="Ensembl.gene",keepMultipleMappings=TRUE, duplicateRemoverMethod="max",orderAbsValue=FALSE)
    plot.object <- data_frame(score=gseaScores(gsca@geneList, gsca@listOfGeneSetCollections[[1]][[1]], exponent=1, mode="graph")$runningScore,position=gseaScores(gsca@geneList, gsca@listOfGeneSetCollections[[1]][[1]], exponent=1, mode="graph")$positions)
    plot.object$x <- 1:nrow(plot.object)
    plot.object
}

df.0 <- data.frame(window=rev(c("0_144","24_72","24_96","72_144","96_144")),xmin=rev(c(0,24,24,72,96)),xmax=rev(c(144,72,96,144,144)),ymin=(c(0,0.5,1,1.5,2)),ymax=(c(0,0.5,1,1.5,2)+0.45))
df.0 <- mutate(df.0,window=factor(window,levels=c("0_144","24_72","24_96","72_144","96_144"),ordered=TRUE))

ggplot(df.0,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=window))+
    geom_rect(alpha=0.7)+
    scale_fill_brewer(palette="Set1")+
    mytheme.model.fit.paper+
    scale_y_continuous("",breaks=c(0,0.5,1,1.5,2)+0.225,labels=rev(c("0-144","24-72","24-96","72-144","96-144")))+
    scale_x_continuous("Time [h]",breaks=c(0,24,72,96,144),labels=c("0h","24h","72h","96h","144h"))+
    guides(fill=FALSE)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle("Window treatment")+
    theme(plot.title=element_text(size=10))+
    geom_vline(xintercept=c(24,72,96),linetype="dotted")+
    theme(plot.background=element_blank())+
    theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+
    theme(axis.title.x=element_blank())->
    p.window.schematic


window.data.96 <- read_csv(file=vpa.acute.study.data)
window.data.96 %>%
    filter(t==144) %>%
    group_by(gene) %>%
    mutate(fc=mny-mny[conc==0]) %>%
    filter(!is.na(start) & !is.na(end)) ->
    window.data.96

window.data.rest <- read_csv(file=vpa.window.data)
window.data.rest %>%
    filter(t==144) %>%
    group_by(gene) %>%
    mutate(fc=mny-mny[conc==0]) %>%
    filter(!is.na(start) & !is.na(end)) ->
    window.data.rest


window.data.all <- bind_rows(window.data.rest,window.data.96)

window.data.all %>%
    ungroup %>%
    filter(gene %in% regulated.genes) %>%
    unite("sample",start,end) %>%
    select(sample,gene,fc) %>%
    rename(value=fc) ->
    gsea.input.tibble

gsea.input.tibble %>%
    filter(gene %in% all.genes) %>%
    group_by(sample) %>%
    do(ll=deframe(.[c("gene","value")])) %$%
    setNames(ll,sample) ->
    gsea.list.in

gmt.object <- getGmt(neural.crest.gmt.file)
gene.set.list <- geneIds(gmt.object)


lapply(names(gsea.list.in),provide.gsea.plot) %>%
    bind_rows(.id="window") %>%
    mutate(label=names(gsea.list.in)[as.numeric(window)])->
    plot.object

ggplot(plot.object,aes(x=x,y=score,group=window,colour=label))+
    geom_hline(yintercept=0,colour="darkgray")+
    geom_segment(aes(x=x,xend=x,y=0+(5-as.numeric(window))/10-0.8,yend=0.8*position/10+(5-as.numeric(window))/10-0.8,group=window),alpha=0.5,size=0.5)+
    geom_line()+
    mytheme.model.fit.paper+
    scale_colour_brewer(palette="Set1")+
    guides(colour=FALSE)+
    scale_y_continuous("Neural Crest\nEnrichment score",breaks=c(0,0.5))+
    scale_x_continuous("Genes ranked by fold change")->
    p.timewindow

p.timewindow <- p.timewindow+theme(plot.margin = unit(c(5, 60, 5, 5), "pt"))


data.pval <- read_csv(gsea.out.nc)

data.pval %>%
    mutate(L1=factor(L1,levels=rev(c("0_144","24_72","24_96","72_144","96_144")),ordered=TRUE)) %>%
    mutate(pval.lab=ifelse(pvals==0,paste("p < ",sprintf("%03.4f",1/data.pval$no.reps[1]),sep=""),paste("p = ",format(pvals,digits=2)))) %>%
    ggplot(aes(x=1,y=L1,label=pval.lab))+
    geom_text(hjust=0,size=3)+
    mytheme.model.fit.paper+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.line.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.line.y=element_blank())+
    theme(panel.background=element_blank(),plot.background=element_blank())->
    p.pval.annot


######################################################################################################
########### Putting everything together
###########
######################################################################################################

row.1 <- plot_grid(p.gsea.all+theme(plot.margin = unit(c(10,5,5,5), "pt")),p.gsea.dev+theme(plot.margin = unit(c(10,5,5,5), "pt")),p.wnt.heatmap+theme(plot.margin = unit(c(10,5,5,5), "pt")),nrow=1,rel_widths=c(0.4,0.6,0.35),labels=c("A","B","C"))
row.2 <- plot_grid(p.kinetic.pathways+theme(plot.margin = unit(c(5, 10, 5, 5), "pt")),p.dose.pathways+theme(plot.margin = unit(c(5, 10, 5, 5), "pt")),p.timewindow,nrow=1,rel_widths=c(0.75,0.75,1),labels=c("D","E","F"))
row.2 <- row.2+draw_plot(p.window.schematic, 0.825, 0.63, 0.15, 0.36)+draw_plot(p.pval.annot, 0.84, 0.15, 0.15, 0.31)

pp <- plot_grid(row.1,row.2,nrow=2)

save_plot(filename=out.file,plot=pp,base_height=5.46,base_aspect_ratio=1.8) #

