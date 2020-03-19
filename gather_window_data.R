library(reshape2)
library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
gene.list.file=args[1]
vpa.window.96 <- args[2]
vpa.window.rest <- args[3]
output.file=args[4]
no.top.regulated=as.numeric(args[5])

gene.list=read_csv(file=gene.list.file,col_names=FALSE) %$% X1


window.data.acute <- read_csv(file=vpa.window.96)
window.data.acute %>%
    filter(t==144) %>%
    group_by(gene) %>%
    mutate(fc=mny-mny[conc==0]) %>%
    filter(!is.na(start) & !is.na(end)) ->
    window.data.acute

window.data <- read_csv(file=vpa.window.rest)
window.data %>%
    filter(t==144) %>%
    group_by(gene) %>%
    mutate(fc=mny-mny[conc==0]) %>%
    filter(!is.na(start) & !is.na(end)) ->
    window.data


window.data.all <- bind_rows(window.data,window.data.acute)

window.data.all %>%
    filter(gene %in% gene.list[1:no.top.regulated]) %>%
    select(start,end,gene,fc) %>%
    unite("start_end",start,end) %>%
    rename(sample=start_end,gene=gene,value=fc) %>%
    write_csv(path=output.file)
