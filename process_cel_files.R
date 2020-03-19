replace.filenames <- function(filenames){
    filenames <- str_replace_all(filenames, "CEL", "")
    filenames <- str_replace_all(filenames, "ND", "")
    filenames <- str_replace_all(filenames, "æ", "")
    filenames <- str_replace_all(filenames, "µ", "")
   filenames <- trimws(filenames)
    filenames <- gsub("[^[:alnum:] ]", "", filenames)
    filenames <- str_replace_all(filenames, "Plus22", "Plus2")
    filenames <- str_replace_all(filenames, "D120878", "D1208")
    filenames <- str_replace_all(filenames, "D121079", "D1210")
    filenames <- str_replace_all(filenames, "D121180", "D1211")
    filenames <- str_replace_all(filenames, "D121481", "D1214")
    filenames <- str_replace_all(filenames, "UKN1102sattimeC2untrD1227102", "UKN1102sattimeC2untrD1227")
    filenames <- str_replace_all(filenames, "UKN1103sattimeC2untrD1228103", "UKN1103sattimeC2untrD1228")
    filenames <- str_replace_all(filenames, "UKN1104sattimeC2untrD1229104", "UKN1104sattimeC2untrD1229")
    filenames <- str_replace_all(filenames, "UKN1105sattimeC2untrD12TW3105", "UKN1105sattimeC2untrD12TW3")
    filenames <- str_replace_all(filenames, "UKN1114sattimeC3untrD1227114", "UKN1114sattimeC3untrD1227") 
    filenames <- str_replace_all(filenames, "UKN1115sattimeC3untrD1228115", "UKN1115sattimeC3untrD1228") 
    filenames <- str_replace_all(filenames, "UKN1116sattimeC3untrD1229116", "UKN1116sattimeC3untrD1229") 
    filenames <- str_replace_all(filenames, "UKN1117sattimeC3untrD12TW3117", "UKN1117sattimeC3untrD12TW3")
    filenames <- str_replace_all(filenames, "UKN1126sattimeC4untrD1227126", "UKN1126sattimeC4untrD1227") 
    filenames <- str_replace_all(filenames, "UKN1127sattimeC4untrD1228127", "UKN1127sattimeC4untrD1228") 
    filenames <- str_replace_all(filenames, "UKN1128sattimeC4untrD1229128", "UKN1128sattimeC4untrD1229") 
    filenames <- str_replace_all(filenames, "UKN1129sattimeC4untrD12TW3129", "UKN1129sattimeC4untrD12TW3")
    }



translate.features <- function(x){
  lbr.string <- paste(annotation(x),".db",sep="")
  library(lbr.string,character.only=TRUE)
  eval(parse(text=paste("map.obj=",paste(annotation(x),"ENSEMBL",sep=""),sep="")))
  print(map.obj)
  mapped_genes <- mappedkeys(map.obj)
                                        # Convert to a list
  xx <- as.list(map.obj[mapped_genes])
  rn <- rownames(exprs(x))[rownames(exprs(x)) %in% mapped_genes]
  expression.vals <- exprs(x)[rn,]
  rownames(expression.vals) <- xx[rn]
  expression.vals
}


compute.eset <- function(all.cels){
    master.batch <- oligo::read.celfiles(filenames=all.cels)
    eset <- oligo::rma(master.batch)
    annotation(eset) <- "hgu133plus2"
    eset <- translate.features(eset)
    rntb <- table(sapply(rownames(eset),nchar))
    normal.length.rowname <- as.numeric(names(rntb[which.max(rntb)]))
    eset=eset[sapply(rownames(eset),nchar)==normal.length.rowname,]
    eset <- aggregate(eset,by=list(rownames(eset)),mean,na.rm=TRUE)
    rownames(eset) <- eset[,"Group.1"]
    eset
}


library(affy)
library(oligo)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
input.dir=args[1]
annotation.file=args[2]
output.file=args[3]

all.files <- do.call("c",lapply(input.dir,list.files,pattern="*.CEL",full.names=T))
all.names <- do.call("c",lapply(input.dir,list.files,pattern="*.CEL",full.names=F))

sysdt.annot <- read_csv(file=annotation.file)


all.names <- replace.filenames(all.names)
sysdt.annot$"File (new name)" <- replace.filenames(sysdt.annot$"File (new name)")
colnames(sysdt.annot)[colnames(sysdt.annot)=="File (new name)"] <- "sample"

tibble(file=all.files,sample=all.names) %>%
    left_join(sysdt.annot) %>%
    filter(!is.na(`time of differentiation [h]`)) ->
    sample.translation.table


eset <- compute.eset(all.files)

colnames(eset) <- replace.filenames(colnames(eset))

eset %>%
    as_tibble %>%
    dplyr::rename(gene=Group1) %>%
    gather("sample","value",-1) %>%
    left_join(sample.translation.table) %>%
    dplyr::rename(t=`time of differentiation [h]`,protocol=`differentiation protocol`,cell_line=`cell line`,treatment_start=`treatment_start [h]`,treatment_end=`treatment_end [h]`,concentration=`compound_concentration [M]`) %>%
    dplyr::select(gene,sample,value,t,replicate,protocol,cell_line,compound,study,treatment_start,treatment_end,concentration) ->
    eset.annotated

write_csv(eset.annotated,path=output.file)
