# oPool<sup>+</sup> display results and experimental data

## Introduction
This folder contains the scripts for oPool<sup>+</sup> display screen result/experimental data processing, analysis, and plotting.

## Folders
* [./experimental_data](./experimental_data): Experimental data of oPool<sup>+</sup> display validation and functional charaterizations
* [./graph](./graph): All plots generated in this study
* [./oPool_results](./oPool_results): Processed screening data and final results
* [./ref_files](./ref_files): Reference files used in this study
* [./script](./script): All custom scripts used in this study

# Step-by-step breakdown of oPool<sup>+</sup> display analysis

## Note
All scripts were executed at oPool_display/oPool_display/ level.

## Input files
* [./ref_files/300lib_Abs.csv](./ref_files/300lib_Abs.csv): Table S1, information of selected antibodies
* [./ref_files/lib_ref.tsv](./ref_files/300lib.tsv): Reference sequences of the natively paired antibody design
* [./ref_files/neg_abs_list.tsv](./ref_files/neg_abs_list.tsv): List of the 30 HA head antibodies (negative controls)
* [./ref_files/202412_sample_name.tsv](./ref_files/202412_sample_name.tsv): Sample names of PacBio sequencing files
* Raw read (PacBio CCS) files in fastq format from NIH SRA database [BioProject PRJNA1150188](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1150188)
* [./experimental_data/](./experimental_data/): Experimental data (validation, structural & functional characterization)

## Analysis of mRNA display results

### Data preparation
(Before processing, adjust the file name of H1 and H3 stem screens according to the sample names table)
1. Filter the CCS reads based on quality and ROI, then trim the adaptors   
``python3 script/filter_from_fastq.py``
    - Input files: PacBio CCS reads
    - Output files will be placed in a folder named fastq_filtered/

### Counting reads, identify scFvs, and calculate enrichment for both replicates
1. Counting unique reads, rename samples
``python3 script/oPool_analysis/processing/fastq2count.py``   
    - Input files:
      - Merged read files in fastq_merged/ folder
      - [./ref_files/202412_sample_name.tsv](./ref_files/202412_sample_name.tsv)   
    - Output files:
      - [./oPool_result/processing/202412_mut_nuc_count.tsv](./oPool_result/processing/202412_mut_nuc_count)

``python3 script/oPool_analysis/processing/split_count_df.py``
    - Input files:
      - [./oPool_result/processing/202412_mut_nuc_count.tsv](./oPool_result/processing/202412_mut_nuc_count)
    - Output files:
      - split count files for different experiments

2. Identify natively paired scFvs with no mutations, then calculate frequncy and enrichment
``python3 script/oPool_analysis/processing/identify_scfv_300lib.py``   
    - Input files:
      - [./oPool_result/nuc_count_files/HA_stem.tsv](./oPool_result/nuc_count_files/HA_stem.tsv) 
      - [./oPool_result/nuc_count_files/Full_HA.tsv](./oPool_result/nuc_count_files/Full_HA.tsv)       
      - [./oPool_result/nuc_count_files/Full_HA_CR9114_competition.tsv](./oPool_result/nuc_count_files/Full_HA_CR9114_competition.tsv)
      - [./ref_files/lib_ref.csv](./ref_files/lib_ref.csv)
    - Output files: 
      - [./oPool_result/enrichment/HA_stem_enrich.tsv](./oPool_result/enrichment/HA_stem_enrich.tsv) 
      - [./oPool_result/enrichment/Full_HA_enrich.tsv](./oPool_result/enrichment/Full_HA_enrich.tsv)       
      - [./oPool_result/enrichment/Full_HA_CR9114_competition_enrich.tsv](./oPool_result/enrichment/Full_HA_CR9114_competition.tsv)    

### Assembly assessment
1. Plot scFv frequency correlation between two assembly replicates
``Rscript script/plot_QC_input.R``   
    - Input files:
      - [./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv](./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv)
    - Output files:
      - [./graph/input_QC.png](./graph/input_QC.png)

### Screening result assessment
1. Plot scFv enrichment correlation between two replicates for each screen
``Rscript script/oPool_analysis/QC/plot_QC_H1_screen.R`` 
``Rscript script/oPool_analysis/QC/plot_QC_H3_screen.R`` 
``Rscript script/oPool_analysis/QC/plot_QC_H1_H3.R`` 
``Rscript script/oPool_analysis/QC/plot_QC_Full_HAs.R`` 
``Rscript script/oPool_analysis/QC/plot_QC_Full_HA_wCompetition.R``
    - Input files:
      - [./oPool_result/enrichment/HA_stem_enrich.tsv](./oPool_result/enrichment/HA_stem_enrich.tsv) 
      - [./oPool_result/enrichment/Full_HA_enrich.tsv](./oPool_result/enrichment/Full_HA_enrich.tsv)       
      - [./oPool_result/enrichment/Full_HA_CR9114_competition_enrich.tsv](./oPool_result/enrichment/Full_HA_CR9114_competition.tsv)   
    - Output files:
      - Correlation plots in graph/oPool_analysis/QC     

### Screening results
1. Plot binding score heatmap
``Rscript script/oPool_analysis/result_analysis/plot_heatmap.R``
    - Input files:
      - [./oPool_result/enrichment/Full_HA_enrich.tsv](./oPool_result/enrichment/Full_HA_enrich.tsv)
      - [./oPool_result/enrichment/Full_HA_CR9114_competition_enrich.tsv](./oPool_result/enrichment/Full_HA_CR9114_competition_enrich.tsv)  
      - [./oPool_result/enrichment/HA_stem_enrich.tsv](./oPool_result/enrichment/HA_stem_enrich.tsv)
    - Output files:
       - [./oPool_result/enrichment/table_s3_1.tsv](./oPool_result/enrichment/table_s3_1.tsv)     
       - [./oPool_result/enrichment/table_s3_1.tsv](./oPool_result/enrichment/table_s3_2.tsv)   
       - [./oPool_result/enrichment/table_s3_1.tsv](./oPool_result/enrichment/table_s3_3.tsv)   
       - [./oPool_result/enrichment/table_s3_1.tsv](./oPool_result/enrichment/table_s3_4.tsv)  
       - [./graph/oPool_heatmap.png](./graph/oPool_heatmap.png)  

2. Cutoff based filtering
``Rscript script/oPool_analysis/result_analysis/cutoff_based_filtering.R``
    - Input files:
      - [./oPool_result/enrichment/combined_enrichment.tsv](./oPool_result/enrichment/combined_enrichment.tsv)
    - Output files:
      - filtered antibody hits for each antigen in oPool_result/filtered_hits

3. Compile the final results table (Table S4)
``python3 script/oPool_analysis/result_analysis/antibody_profiling.R``       
    - Input files:
      - [./oPool_result/enrichment/combined_enrichment.tsv](./oPool_result/enrichment/combined_enrichment.tsv)
    - Output files:
      - [./oPool_result/filtered_hits/unique_antibodies_info.tsv](./oPool_result/filtered_hits/unique_antibodies_info.tsv)

4. Plot competition indices
``Rscript script/oPool_analysis/result_analysis/plot_competition_index.R``
    - Input files:
      - filtered antibody hits for each antigen in oPool_result/filtered_hits
    - Output files:
      - competition index bar plots in graph/oPool_analysis/competition
      - [./oPool_result/enrichment/competition_index.tsv](./oPool_result/enrichment/competition_index.tsv)

## Analysis of experimental validation results

### Hits validation via BLI
1. Plot quantitation data
``Rscript script/plot_scFv_quantitation_standard_curve.R`   
    - Input files:
      - [./result/BLI_data/scFv_quantitation.tsv](./result/BLI_data/scFv_quantitation.tsv)
    - Output files:
      - [./graph/scFv_quantitation_curve.png](./graph/scFv_quantitation_curve.png)

2. Compile raw data for plotting
``python3 script/compile_BLI_result_batch.py`   
    - Input files:
      - [./result/BLI_data/](./result/BLI_data/):raw BLI data
    - Output files:
      - [./result/BLI_compile/](./result/BLI_compile/): all compiled BLI data

3. Plot kinetics data
``Rscript script/plot_BLI_binding_batch.R`   
    - Input files:
      - [./result/BLI_compile/](./result/BLI_compile/): all compiled BLI data
    - Output files:
      - [./graph/BLI_sensorgram/](./graph/BLI_sensorgram/): all sensorgrams

4. Plot validation heatmap
``Rscript script/plot_validation_heatmap.R`   
    - Output files:
      - [./graph/H1_validation_heatmap.png](./graph/H1_validation_heatmap.png)
      - [./graph/H3_validation_heatmap.png](./graph/H3_validation_heatmap.png)

## Characterization of AG11-2F01 and 16.ND.92

### Structural analyses 
1. Plot structure overviews of AG11-2F01 and 16.ND.92
``pymol script/overview.pml`   
    - Input files:
      - [./data/PDB/SI06HA_2F01.pdb](./data/PDB/SI06HA_2F01.pdb)
      - [./data/PDB/SI06HA_16ND92.pdb](./data/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/PDB/2F01_overview.png](./graph/PDB/2F01_overview.png)
      - [./graph/PDB/16ND_overview.png](./graph/PDB/16ND_overview.png)

2. Plot epitopes of AG11-2F01 and 16.ND.92
``pymol script/epitope.pml`   
    - Input files:
      - [./data/PDB/SI06HA_2F01.pdb](./data/PDB/SI06HA_2F01.pdb)
      - [./data/PDB/SI06HA_16ND92.pdb](./data/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/PDB/2F01_epitope.png](./graph/PDB/2F01_epitope.png)
      - [./graph/PDB/16ND_epitope.png](./graph/PDB/16ND_epitope.png)

3. Plot CDRH3 interactions of AG11-2F01 and 16.ND.92
``pymol script/interact_CDRH3.pml`   
    - Input files:
      - [./data/PDB/SI06HA_2F01.pdb](./data/PDB/SI06HA_2F01.pdb)
      - [./data/PDB/SI06HA_16ND92.pdb](./data/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/PDB/2F01_CDRH3.png](./graph/PDB/2F01_CDRH3.png)
      - [./graph/PDB/16ND_CDRH3.png](./graph/PDB/16ND_CDRH3.png)

4. Plot light chain interactions  of AG11-2F01 and 16.ND.92
``pymol script/interact_LC.pml`   
    - Input files:
      - [./data/PDB/SI06HA_2F01.pdb](./data/PDB/SI06HA_2F01.pdb)
      - [./data/PDB/SI06HA_16ND92.pdb](./data/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/PDB/2F01_LC.png](./graph/PDB/2F01_LC.png)
      - [./graph/PDB/16ND_LC.png](./graph/PDB/16ND_LC.png)

5. Plot BSA of IGHD3-3 HA stem antibodies
``Rscript script/plot_BSA_bar_chart.R`   
    - Input files:
      - [./data/bsa_percentage.tsv](./data/bsa_percentage.tsv)
    - Output files:
      - [./graph/BSA_stacked_bar_chart.png](./graph/BSA_stacked_bar_chart.png)

6. Plot IGHD3-3 contribution to VH paratopes
``Rscript script/plot_3-3_BSA.R`   
    - Input files:
      - [./data/3-3_BSA_percentage.tsv](./data/3-3_BSA_percentage.tsv)
    - Output files:
      - [./graph/BSA_IGHD3-3_bar_plot.png](./graph/BSA_IGHD3-3_bar_plot.png)

7. Plot CDRH3 overlays of IGHD3-3 HA stem antibodies
``pymol script/interact_LC.pml`   
    - Input files:
      - [./data/PDB/SI06HA_2F01.pdb](./data/PDB/SI06HA_2F01.pdb)
      - [./data/PDB/SI06HA_16ND92.pdb](./data/PDB/SI06HA_16ND92.pdb)
      - [./data/PDB/1G05_mono.pdb](./data/PDB/1G05_mono.pdb)
      - [./data/PDB/39.29_mono.pdb](./data/PDB/39.29_mono.pdb)
      - [./data/PDB/56.a.09_mono.pdb](./data/PDB/56.a.09_mono.pdb)
      - [./data/PDB/429_B01_mono.pdb](./data/PDB/429_B01_mono.pdb)
      - [./data/PDB/MEDI8852_mono.pdb](./data/PDB/MEDI8852_mono.pdb)
      - [./data/PDB/SIA28_mono.pdb](./data/PDB/SIA28_mono.pdb)
    - Output files:
      - [./graph/PDB/CDRH3_overlay.png](./graph/PDB/CDRH3_overlay.png)

### Functional characterizations
8. Plot ELISA result heatmap
``Rscript script/plot_EC50.R`   
    - Output files:
      - [./graph/ELISA_EC50_heatmap.png](./graph/ELISA_EC50_heatmap.png)

9. Plot micro-neutralization result heatmap
``Rscript script/plot_IC50.R`   
    - Output files:
      - [./graph/IC50_heatmap.png](./graph/IC50_heatmap.png)

10. Plot in vivo experiment data
``Rscript script/plot_IC50.R` 
    - Input files:
      - [./data/invivo_weight_loss.tsv](./data/invivo_weight_loss.tsv)
      - [./data/invivo_survival.tsv](./data/invivo_survival.tsv)
      - [./data/invivo_lung_titer.tsv](./data/invivo_lung_titer.tsv)
    - Output files:
      - [./graph/invivo_weight_loss_2F01.png](./graph/invivo_weight_loss_2F01.png)
      - [./graph/invivo_weight_loss_16ND92.png](./graph/invivo_weight_loss_16ND92.png)
      - [./graph/invivo_survival_2F01.png](./graph/invivo_survival_2F01.png)
      - [./graph/invivo_survival_16ND92.png](./graph/invivo_survival_16ND92.png)
      - [./graph/invitro_lung_titer_Vero_2F01.png](./graph/invitro_lung_titer_Vero_2F01.png)
      - [./graph/invitro_lung_titer_Vero_16ND92.png](./graph/invitro_lung_titer_Vero_16ND92.png)


