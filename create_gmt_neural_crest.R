library(GSEABase)
library(org.Hs.eg.db)
library(GO.db)
library(EnrichmentBrowser)

args = commandArgs(trailingOnly=TRUE)
input.file=args[1]
output.file=args[2]

full.msigdb <- input.file
foo <- getGmt(full.msigdb)
broad.sig.list <- geneIds(foo)

gs <- GeneSetCollection(GeneSet(broad.sig.list[["LEE_NEURAL_CREST_STEM_CELL_UP"]],setName="Neural_Crest"))
Tfile <- file(output.file, "w+")
toGmt(gs,con=Tfile)
close(Tfile)
