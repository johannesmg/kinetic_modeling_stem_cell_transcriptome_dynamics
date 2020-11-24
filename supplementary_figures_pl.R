args = commandArgs(trailingOnly=TRUE)
fit.file=args[1]
all.gene.file=args[2]
symbol.conversion.file=args[3]
n.top.regulated=as.numeric(args[4])
pl.directory=args[5]
no.cores=args[6]
loess.span=as.numeric(args[7])
parameter.boundary.file=args[8]

source("common_setup.R")


model.based.peak.time=function(x){
    x %>%
        ungroup %>%
        group_by(message,run,dt) %>%
        do(tibble(t=seq(0,144,1),yutr=bateman.model(seq(0,144,1),conc=rep(0,length(seq(0,144,1))),parms=as.numeric(.[1,par.names])),ytr=bateman.model(seq(0,144,1),conc=rep(0.6,length(seq(0,144,1))),parms=as.numeric(.[1,par.names])))) %>%
        gather("type","value",yutr:ytr) %>%
        group_by(message,run,type) %>%
        summarise(peak.time=t[which.max(value)],max.exp=max(value),min.exp=min(value)) ->
        peak.time.out
    peak.time.out
}

######################################################################################################
########### Load expression data and determine data based peak times
###########
######################################################################################################


fit.data=read_csv(fit.file)

par.names <- c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")


fit.data %>%
    filter(conc == 0) %>%
    group_by(gene) %>%
    do(ll=data_frame(t=seq(0,144,0.1),y=predict(loess(mny~t,data=.,span=loess.span),newdata=seq(0,144,0.1)))) %$%
    setNames(ll,gene) ->
    interpolated.data


gene.maximal.reg <- enframe(unlist(lapply(interpolated.data,function(x){x$t[which.max(x$y)]}))) %>%
    rename(gene=name,peaktime=value)

peaking.intermediate=gene.maximal.reg[gene.maximal.reg$peaktime>0 & gene.maximal.reg$peaktime< 144,]$gene

######################################################################################################
########### Load PL data
###########
######################################################################################################


options(readr.num_columns = 0)

in.files=list.files(paste(pl.directory,"/all",sep=""),full.names=T,pattern="checked*")
optimality=bind_rows(mclapply(in.files,read_csv,mc.cores=no.cores))

optimality %>%
    select(gene,status) %>%
    mutate(filename=paste("pl_outgene_",gene,".csv",sep="")) %>%
    mutate(path=ifelse(status=="rerun_not_necessary","all/","second/")) ->
    optimality

file.names.to.load=lapply(1:nrow(optimality),function(x){paste(pl.directory,optimality[x,]$path,optimality[x,]$filename,sep="")})


full.pl=mclapply(file.names.to.load,read_csv,guess_max=4e3,mc.cores=no.cores)
full.pl=bind_rows(full.pl)

######################################################################################################
########### Confidence intervals for gon,goff,n,dt
###########
######################################################################################################


full.pl[,c("gene","method","message",par.names,"parameter","parameter.value")] %>%
    filter(gene %in% regulated.genes) ->
    full.pl.not.spread

full.pl.not.spread %>%
    mutate(parameter=ifelse(parameter=="g1","gon",parameter)) %>%
    mutate(parameter=ifelse(parameter=="g2","goff",parameter)) %>%
    group_by(gene,method,parameter) %>%
    summarise(mn=min(parameter.value),mx=max(parameter.value),opt=parameter.value[grepl("pre pl",message)]) %>%
    mutate(rel.mn=mn/opt,rel.mx=mx/opt) ->
    pl.borders

pl.borders %>%
    mutate(range=mx/mn) %>%
    group_by(gene,parameter) %>%
    summarise(range=max(range)) %>%
    filter(parameter %in% c("gon","goff","n","dt")) %>%
    ungroup %>%
    left_join(gene.maximal.reg) %>%
    filter(peaktime>0 & peaktime< 144) %>%
    mutate(peaktime=cut(peaktime,seq(0,144,12)))  %>%
    group_by(parameter,peaktime) %>%
    do(y=enframe(quantile(log10(.$range),c(0.25,0.5,0.75),na.rm=T))) %>%
    unnest %>%
    spread(name,value) %>%
    ggplot(aes(x=peaktime,y=`50%`,ymin=`25%`,ymax=`75%`)) +
    geom_crossbar()+
    facet_wrap(~parameter,scales="free")+
    theme_minimal()+
    theme(axis.text.x=element_text(angle=90,vjust=0.5))+
    scale_x_discrete("Gene peak time (data based) [h]")+
    scale_y_continuous("Log10 (right CI border / left CI border)") ->
    p1

save_plot(p1,filename = paste("./figures/supplement_pl_g1_g2_vs_peak_time",as.character(n.top.regulated),".pdf",sep=""),base_height=5)

rm(full.pl.not.spread)

######################################################################################################
########### Uncertainty of approximate analytically predicted peak time
###########
######################################################################################################


full.pl[,c("gene","method","message",par.names,"parameter","parameter.value")] %>%
    filter(gene %in% regulated.genes) %>%
    group_by(gene,method) %>%
    mutate(run=1:n()) %>%
    ungroup %>%
    gather("other.par","value",g1p:dt) %>%
    mutate(value=ifelse(parameter==other.par,parameter.value,value)) %>%
    select(gene,method,message,run,other.par,value) %>%
    spread(other.par,value) ->
    full.pl



full.pl %>%
    mutate(peak.time=(n-1)/g1+1/g2-dt) %>%
    group_by(gene,method) %>%
    summarise(peak.time.max=max(peak.time),peak.time.min=min(peak.time))->
    par.out


gene.maximal.reg %>%
    right_join(par.out) %>%
    mutate(df=(peak.time.max-peak.time.min)) %>%
    group_by(peaktime,gene) %>%
    summarise(df=max(df)) %>%
    ungroup %>%
    filter(peaktime>0 & peaktime<=108) %>%
    mutate(peaktime=cut(peaktime,seq(0,144,12),include.lowest=T)) %>%
    group_by(peaktime) %>%
    summarise(md=median(df)) %>%
    ggplot(aes(x=peaktime,y=md))+
    geom_bar(stat="identity")+
    theme_minimal()+
    scale_x_discrete("Gene peak time (data based) [h]")+
    scale_y_continuous("Median (n-1)/g_on+1/g_off-dt 95% CI")->
    p1

save_plot(p1,filename=paste("./figures/supplement_pl_peaktime_parameters_ci_",as.character(n.top.regulated),".pdf",sep=""),base_height=5)

######################################################################################################
########### Compute amplitude and peak time from PL data
###########
######################################################################################################



full.pl %>%
    unite("gene_method",gene,method) %>%
    group_by(gene_method) ->
    full.pl.df.in

full.pl.df.in %>%
    group_split() %>%
    set_names(unlist(group_keys(full.pl.df.in)))  ->
    full.pl.list

amplitude.peaktime.out=mclapply(full.pl.list[],model.based.peak.time,mc.cores=no.cores)


bind_rows(amplitude.peaktime.out,.id="gene_method") %>%
    separate(gene_method,into=c("gene","method"),sep="_") %>%
    select(gene,method,message,run,type,peak.time,max.exp,min.exp) ->
    amplitude.peaktime.out.frame


######################################################################################################
########### Amplitude change significance
###########
######################################################################################################



amplitude.peaktime.out.frame %>%
    ungroup %>%
    filter(grepl("pre pl",message)) %>%
    distinct(gene,method,type,.keep_all=T) %>%
    mutate(diff.exp=max.exp-min.exp) %>%
    select(gene,method,run,type,diff.exp) %>%
    spread(type,diff.exp) %>%
    mutate(diff.amp=ytr-yutr) ->
    optimal.parameter.set.amplitude.diff
    

amplitude.peaktime.out.frame %>%
    mutate(diff.exp=max.exp-min.exp) %>%
    select(gene,method,run,type,diff.exp) %>%
    spread(type,diff.exp) %>%
    mutate(diff.amp=ytr-yutr) %>%
    group_by(gene) %>%
    summarise(max.amp=max(diff.amp),min.amp=min(diff.amp)) %>%
    mutate(status=ifelse((min.amp >0 & max.amp >0)|(min.amp<0 & max.amp <0),"significant","not significant")) ->
    amplitude.status

optimal.parameter.set.amplitude.diff %>%
    left_join(amplitude.status) %>%
    ggplot(aes(x=diff.amp,group=status))+
    geom_histogram()+
    geom_vline(xintercept=0,colour="red",alpha=0.5,size=1)+
    facet_wrap(~status)+
    scale_x_continuous("Log2 Amplitude shift")+
    theme_minimal()->
    p1


save_plot(p1,filename="./figures/supplement_pl_amplitude_change_ci_corrected.pdf",base_height=5)

######################################################################################################
########### Peak time change significance
###########
######################################################################################################


amplitude.peaktime.out.frame %>%
    group_by(gene,type) %>%
    summarise(ymin=min(peak.time)[1],ymax=max(peak.time)[1],y=peak.time[grepl("pre pl",message)][1]) ->
    peak.info

peak.info %>%
    select(gene,type,y) %>%
    spread(type,y) %>%
    filter(yutr>0 & yutr<144) %>%
    mutate(mean.peak=(ytr+yutr)/2,diff.peak=ytr-yutr) ->
    peak.diff.optimal.parameters

peak.info %>%
    gather("extremum","value",ymin:ymax) %>%
    group_by(gene) %>%
    summarise(max.diff=max(range(value[type=="ytr"]))-min(range(value[type=="yutr"])),min.diff=min(range(value[type=="ytr"]))-max(range(value[type=="yutr"]))) %>%
    mutate(status=ifelse(max.diff>0 & min.diff >0, "right",ifelse(max.diff<0 & min.diff<0,"left","unclear")))->
    status.info


peak.diff.optimal.parameters %>%
    left_join(select(status.info,gene,status)) %>%
    mutate(status=ifelse(status!="unclear","significant","not significant")) %>%
    ggplot(aes(x=mean.peak,y=diff.peak,group=status,colour=status))+
    geom_point()+
    theme_minimal()+
    scale_colour_manual(values=c("darkgray","darkred"))+
    stat_smooth(geom="line",se=F,size=1,alpha=0.5,method="loess",span=0.75)+
    scale_x_continuous("Average peak time")+
    scale_y_continuous("Peak shift upon VPA treatment")->
    p1

save_plot(p1,filename="./figures/supplement_pl_peaktime_change_significance.pdf",base_height=5)


######################################################################################################
########### Example PL for POU4F1
###########
######################################################################################################



pl.gene=read_csv(paste(pl.directory,"/all/pl_outgene_ENSG00000152192.csv",sep=""))

cost.increase=qchisq(0.95,1)
pp=read_csv(parameter.boundary.file)

pp %>%
    gather("parameter","value",g1p:dt) %>%
    mutate(value=log10(value)) ->
    bounds.frame

pl.gene %>%
    group_by(parameter) %>%
    mutate(parameter.value=log10(parameter.value)) %>%
    summarise(mx=max(parameter.value),mn=min(parameter.value)) %>%
    gather("type","observed.value",mx:mn) %>%
    left_join(bounds.frame) %>%
    filter(abs(value-observed.value)<0.05) ->
    bounds.frame

pl.gene %>%
    filter(method=="lbfgsb") %>%
    select(gene,parameter,parameter.value,cost,method,message) %>%
    group_by(gene) %>%
    mutate(cost=cost-min(cost)[1]) %>%
    mutate(parameter.value=log10(parameter.value)) ->
    pl.gene

pl.gene %>%
    ggplot(aes(x=parameter.value,y=cost,group=method))+
    geom_vline(data=bounds.frame,aes(xintercept=value),colour="gray20",linetype="dotted")+
    geom_hline(yintercept = cost.increase,colour="gray20",linetype="dotted")+
    geom_point(size=2,data=filter(pl.gene,message=="pre pl optimum"))+
    geom_line()+
    facet_wrap(~parameter,scales="free_x")+
    theme_minimal()+
    coord_cartesian(ylim=c(0,4.2))+
    scale_x_continuous("Log10 parameter value")+
    scale_y_continuous("Cost increase relative to optimum") ->
    p1

save_plot(plot=p1,"./figures/supplement_pl_example_gene.pdf",base_height = 5)
