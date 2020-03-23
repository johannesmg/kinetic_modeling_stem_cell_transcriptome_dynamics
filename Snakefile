configfile: 'config.yml'

ESET_FILE=config['eset_file']
SYMBOL_CONVERSION_FILE=config['symbol_conversion_file']
PARAMETER_FILE=config['parameter_file']
FULL_MSIGDB=config['full_msigdb']
KEGG_MSIGDB=config['kegg_msigdb']
RAW_DATA_PCR=config['raw_data_pcr']
CEL_DIR=config['cel_dir']
CEL_ANNOTATION_FILE=config['annotation_file']

N_TOP_REGULATED=config['n_top_regulated']
MINSETSIZE_MSIGDB=config['minsetsize_msigdb']
MINSETSIZE_KEGG=config['minsetsize_kegg']
MINSETSIZE_NC=config['minsetsize_nc']
NO_PERMUTATIONS=config['nopermutations']
NO_CORES_GSEA=config['nocoresgsea']
SIG_THRESH_GSEA_KEGG=config['sig_thresh_gsea_kegg']
DATA_TIME_ROW=config['data_time_row']
DATA_CONC_ROW=config['data_conc_row']

TEST_TIME= [x*4 for x in range(0, 36)]
TEST_CONC= [x/40 for x in range(1, 24)]

# A quick test run can be made with this configuration
# TEST_TIME= [x*4 for x in range(0, 2)]
# TEST_CONC= [x/40 for x in range(1, 2)]




rule all:
 input:
   [expand('figures/panel_01_{n_top_regulated}.pdf',n_top_regulated=N_TOP_REGULATED),expand('figures/panel_02_{n_top_regulated}.pdf',n_top_regulated=N_TOP_REGULATED),expand('figures/panel_03_{n_top_regulated}.pdf',n_top_regulated=N_TOP_REGULATED),expand('figures/panel_04_{sig_thresh_gsea_kegg}_{minsetsize_msigdb}_{minsetsize_kegg}_{n_top_regulated}.pdf',n_top_regulated=N_TOP_REGULATED,minsetsize_kegg=MINSETSIZE_KEGG,minsetsize_msigdb=MINSETSIZE_MSIGDB,sig_thresh_gsea_kegg=SIG_THRESH_GSEA_KEGG)]

rule create_expression_set:
  input:
    [CEL_ANNOTATION_FILE]
  output:
    [ESET_FILE]
  shell:
    'Rscript process_cel_files.R {CEL_DIR} {input} {ESET_FILE}'


rule fit_input_data:
  input:
    [ESET_FILE]
  output:
    ['output_data/fit_data_in.csv','output_data/all_changing_genes.csv','output_data/vpa_acute_study_window_treatment_data.csv','output_data/vpa_window_treatment_data.csv']
  shell:
    'Rscript eset_to_fit_input.R {ESET_FILE} ./output_data/fit_data_in.csv'


rule mica_dataset:
  input:
    [SYMBOL_CONVERSION_FILE]
  output:
    'output_data/mica_gene_sets.RData'
  shell:
    'Rscript provide_mica_dataset.R {input} {output}'


rule create_gsea_model_data_input:
  input:
    [PARAMETER_FILE,'output_data/fit_data_in.csv','output_data/all_changing_genes.csv']
  output:
    'output_data/gsea_model_data_in_{n_top_regulated}.csv'
  shell:
    'Rscript create_gsea_input.R {input[0]} {output} {input[1]} {input[2]} {wildcards.n_top_regulated}'

    
rule create_nc_gmt:
  input:
    [FULL_MSIGDB]
  output:
    'output_data/neural_crest.gmt'
  shell:
    'Rscript create_gmt_neural_crest.R {input} {output}'


rule create_window_data:
  input:
    ['output_data/all_changing_genes.csv','output_data/vpa_acute_study_window_treatment_data.csv','output_data/vpa_window_treatment_data.csv']
  output:
    'output_data/window_data_gsea_in_{n_top_regulated}.csv'
  shell:
    'Rscript gather_window_data.R {input[0]} {input[1]} {input[2]} {output} {wildcards.n_top_regulated}'


rule gsea_msigdb:
  input:
    ['output_data/gsea_model_data_in_{n_top_regulated}.csv',FULL_MSIGDB]
  output:
    'output_data/gsea_msigdb_{minsetsize}_{n_top_regulated}_{type}_{time}_{conc}.csv'
  shell:
    'Rscript gsea_fc_point.R {input[0]} {input[1]} {wildcards.minsetsize} {NO_PERMUTATIONS} {NO_CORES_GSEA} {wildcards.type} {wildcards.time} {wildcards.conc} {output}'

rule gsea_kegg:
  input:
    ['output_data/gsea_model_data_in_{n_top_regulated}.csv',KEGG_MSIGDB]
  output:
    'output_data/gsea_kegg_{minsetsize}_{n_top_regulated}_{type}_{time}_{conc}.csv'
  shell:
    'Rscript gsea_fc_point.R {input[0]} {input[1]} {wildcards.minsetsize} {NO_PERMUTATIONS} {NO_CORES_GSEA} {wildcards.type} {wildcards.time} {wildcards.conc} {output}'


rule combinegsea:
    input:
      [expand("output_data/gsea_{{set}}_{{minsetsize}}_{{n_top_regulated}}_model_{time}_{conc}.csv",time=TEST_TIME,conc=TEST_CONC),expand("output_data/gsea_{{set}}_{{minsetsize}}_{{n_top_regulated}}_data_{time}_{conc}.csv",time=DATA_TIME_ROW,conc=0.6),expand("output_data/gsea_{{set}}_{{minsetsize}}_{{n_top_regulated}}_data_{time}_{conc}.csv",time=144,conc=DATA_CONC_ROW)]
    output:
      'output_data/final_gsea_{set}_{minsetsize}_{n_top_regulated}.csv'
    shell:
      """
      grep -hv "set.definition.file" ./output_data/gsea_{wildcards.set}_{wildcards.minsetsize}_{wildcards.n_top_regulated}* >> {output}
      rm ./output_data/gsea_{wildcards.set}_{wildcards.minsetsize}_{wildcards.n_top_regulated}*
      """


rule gsea_nc:
  input:
    ['output_data/window_data_gsea_in_{n_top_regulated}.csv','output_data/neural_crest.gmt']
  output:
    'output_data/gsea_neural_crest_window_out_{n_top_regulated}.csv'
  shell:
    'Rscript gsea_fc_multiple.R {input[0]} {input[1]} {MINSETSIZE_NC} {NO_PERMUTATIONS} 1 {output}'


rule panel1:
  input:
    ['output_data/fit_data_in.csv','output_data/all_changing_genes.csv',SYMBOL_CONVERSION_FILE]
  output:
    'figures/panel_01_{n_top_regulated}.pdf'
  shell:
    'Rscript cowplot_introduction_panel.R {input[0]} {input[1]} {input[2]} {wildcards.n_top_regulated} {output}'

rule panel2:
  input:
    ['output_data/fit_data_in.csv','output_data/all_changing_genes.csv',SYMBOL_CONVERSION_FILE,PARAMETER_FILE]
  output:
    'figures/panel_02_{n_top_regulated}.pdf'
  shell:
    'Rscript cowplot_fit_setup_and_evaluation_panel.R {input[0]} {input[1]} {input[2]} {input[3]} {wildcards.n_top_regulated} {output}'

rule panel3:
  input:
    ['output_data/fit_data_in.csv','output_data/all_changing_genes.csv',SYMBOL_CONVERSION_FILE,PARAMETER_FILE,RAW_DATA_PCR]
  output:
    'figures/panel_03_{n_top_regulated}.pdf'
  shell:
    'Rscript cowplot_validation_panel.R {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {wildcards.n_top_regulated} {output}'

rule panel4:
  input:
    ['output_data/fit_data_in.csv','output_data/all_changing_genes.csv',SYMBOL_CONVERSION_FILE,PARAMETER_FILE,FULL_MSIGDB,'output_data/final_gsea_msigdb_{minsetsize_msigdb}_{n_top_regulated}.csv','output_data/final_gsea_kegg_{minsetsize_kegg}_{n_top_regulated}.csv','output_data/vpa_window_treatment_data.csv','output_data/vpa_acute_study_window_treatment_data.csv','output_data/mica_gene_sets.RData',KEGG_MSIGDB,'output_data/gsea_neural_crest_window_out_{n_top_regulated}.csv','output_data/neural_crest.gmt']
  output:
    'figures/panel_04_{sig_thresh_gsea_kegg}_{minsetsize_msigdb}_{minsetsize_kegg}_{n_top_regulated}.pdf'
  shell:
    'Rscript cowplot_functional_panel.R {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} {wildcards.n_top_regulated} {output} {wildcards.sig_thresh_gsea_kegg} {input[12]}'

