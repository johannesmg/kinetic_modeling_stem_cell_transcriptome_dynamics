# kinetic_modeling_stem_cell_transcriptome_dynamics

These scripts generate figures 1-4 for the article "Kinetic modeling of stem cell transcriptome dynamics to identify regulatory modules of normal and disturbed neuroectodermal differentiation". Here, the GPL (>=2) applies.
 
The Snakefile describes the workflow and can be executed using snakemake, downloadable from https://github.com/snakemake/snakemake. In order to execute snakemake, the local CEL file path (see below) has to be entered into config.yml.

In order to create the figure pdfs, the following additional files are required:

1. Two GMT gene set definition files from the MsigDB (./input_data/c2.all.v6.0.entrez.gmt,./input_data/c2.cp.kegg.v6.0.entrez.gmt)
Both files can be downloaded from https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/6.0/msigdb_v6.0_files_to_download_locally.zip

2. Raw microarray CEL files that were uploaded to GEO under the accession number GSE147270.

The file gsea_functions.R contains modified code from the R package HTSAnalyzeR (https://bioconductor.org/packages/release/bioc/html/HTSanalyzeR.html), licensed under the Artistic-2.0.