## Step-by-step breakdown of oPool<sup>+</sup> display analysis

## Contents
- [Introduction](#introduction)
- [Folders](#Folders)
- [Input files](#Input-files)
- [Analysis of oPool<sup>+</sup> display results](#Analysis-of-oPool<sup>+</sup>-display-results)
- [Analysis of experimental validation results](#Analysis-of-experimental-validation-results)
- [Characterization of AG11-2F01 (PDB ID 9DBX) and 16.ND.92 (PDB ID 9CU7)](#Characterization-of-AG11-2F01-and-16.ND.92)

## Introduction
This folder contains the scripts for oPool<sup>+</sup> display screen result/experimental data processing, analysis, and plotting.

## Folders
* [./experimental_data](./experimental_data): Experimental data of oPool<sup>+</sup> display validation and functional characterizations
* [./graph](./graph): All plots generated in this study
* [./oPool_result](./oPool_result): Processed screening data and final results
* [./ref_files](./ref_files): Reference files used in this study
* [./script](./script): All custom scripts used in this study

## Note
All scripts were executed at this level.

## Input files
* [./ref_files/300lib_Abs.tsv](./ref_files/300lib_Abs.tsv): Table S1, information of selected antibodies
* [./ref_files/lib_ref.csv](./ref_files/lib_ref.csv): Reference sequences of the natively paired antibody design
* [./ref_files/neg_abs_list.tsv](./ref_files/neg_abs_list.tsv): List of the 30 HA head antibodies (negative controls)
* [./ref_files/202412_sample_name.tsv](./ref_files/202412_sample_name.tsv): Sample names of PacBio sequencing files
* Raw read (PacBio CCS) files in fastq format from NIH SRA database [BioProject PRJNA1150188](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1150188)
* [./experimental_data](./experimental_data): Experimental data (validation, structural & functional characterization)

## Analysis of oPool<sup>+</sup> display results

### Data preparation
(Before processing, adjust the numerical ID in the sequencing file name of H1 and H3 stem screens according to the sample names table [./ref_files/202412_sample_name.tsv](./ref_files/202412_sample_name.tsv))

1. Filter the CCS reads based on quality and ROI, then trim the adaptors   
``python3 script/oPool_analysis/processing/filter_from_fastq.py``
    - Input files: PacBio CCS reads
    - Output files will be placed in a folder named fastq_filtered/

### Counting reads, identify scFvs, and calculate enrichment for both replicates
1. Counting unique reads, rename samples
``python3 script/oPool_analysis/processing/fastq2count.py``   
    - Input files:
      - Merged read files in fastq_merged/ folder
      - [./ref_files/202412_sample_name.tsv](./ref_files/202412_sample_name.tsv)   
    - Output files:
      - [./oPool_result/processing/202412_mut_nuc_count.tsv](./oPool_result/processing/202412_mut_nuc_count.tsv)

``python3 script/oPool_analysis/processing/split_count_df.py``
    - Input files:
      - [./oPool_result/processing/202412_mut_nuc_count.tsv](./oPool_result/processing/202412_mut_nuc_count.tsv)
    - Output files:
      - split count files for different experiments

2. Identify natively paired scFvs with no mutations, then calculate frequency and enrichment
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
1. Identify natively paired scFvs with no mutations in each assembly, then calculate frequency and enrichment
``python3 script/oPool_analysis/assembly_QC/identify_scfv_small_pool_QC.py``   
    - Input files:
      - [./oPool_result/nuc_count_files/Assembly_*.tsv] read count files for each assembly
    - Output files:
      - [./oPool_result/assembly_QC/freq](./oPool_result/assembly_QC/freq) frequency tables for each assembly

2. Plot scFv frequency correlation between each assembly replicates
``Rscript script/assembly_QC/plot_small_pool_pcr_freq.R`` 
    - Input files:
      - [./oPool_result/assembly_QC/freq](./oPool_result/assembly_QC/freq) frequency tables for each assembly
    - Output files:
      - [./graph/assembly_QC](./graph/assembly_QC) frequency correlation plots for each assembly, Figure 1E

3. Plot scFv frequency correlation between final library replicates
``Rscript script/assembly_QC/plot_QC_input.R``   
    - Input files:
      - [./oPool_result/enrichment/HA_stem_enrich.tsv](./oPool_result/enrichment/HA_stem_enrich.tsv)
    - Output files:
      - [./graph/assembly_QC/input_QC.png](./graph/assembly_QC/input_QC.png)


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
       - [./oPool_result/enrichment/table_s3_2.tsv](./oPool_result/enrichment/table_s3_2.tsv)   
       - [./oPool_result/enrichment/table_s3_3.tsv](./oPool_result/enrichment/table_s3_3.tsv)   
       - [./oPool_result/enrichment/table_s3_4.tsv](./oPool_result/enrichment/table_s3_4.tsv)  
       - [./graph/oPool_heatmap.png](./graph/oPool_heatmap.png): Figure 2

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
      - competition index bar plots in graph/oPool_analysis/competition, Figure 4A
      - [./oPool_result/enrichment/competition_index.tsv](./oPool_result/enrichment/competition_index.tsv)

## Analysis of experimental validation results

### Binding validation via BLI
1. Plot BLI binding validation sensorgrams
``Rscript script/validation/BLI/plot_oPool_binding_validation.R``  
    - Input files:
      - [./experimental_data/validation/BLI/oPool_binding_validation](./experimental_data/validation/BLI/oPool_binding_validation):raw data
      - [./experimental_data/validation/BLI/oPool_binding_validation/oPool_validation_sample_names.xlsx](./experimental_data/validation/BLI/oPool_binding_validation/oPool_validation_sample_names.xlsx):Sample info for each file.
    - Output files:
      - [./graph/validation/BLI/oPool_binding_validation](./graph/validation/BLI/oPool_binding_validation): sensorgrams by antibody

2. Plot BLI binding validation heatmaps
``Rscript script/validation/BLI/plot_oPool_binding_validation_heatmap.R``
    - Input files:
      - [./experimental_data/validation/BLI/oPool_binding_validation](./experimental_data/validation/BLI/oPool_binding_validation):raw data
      - [./experimental_data/validation/BLI/oPool_binding_validation/oPool_validation_sample_names.xlsx](./experimental_data/validation/BLI/oPool_binding_validation/oPool_validation_sample_names.xlsx):Sample info for each file.
    - Output files:
      - [./graph/validation/BLI/oPool_binding_validation/BLI_binding_response_heatmap.png](./graph/validation/BLI/oPool_binding_validation/BLI_binding_response_heatmap.png): Figure 3A

### Binding validation via ELISA
1. Plot ELISA validation heatmap
``Rscript script/oPool_analysis/result_analysis/plot_competition_index.R``
    - Input files:
      - [./experimental_data/validation/ELISA/ELISA_Validation_Results.xlsx](./experimental_data/validation/ELISA/ELISA_Validation_Results.xlsx): OD450 data
    - Output files:
      - [./graph/validation/ELISA/oPool_validation_heatmap_ELISA.png](./graph/validation/ELISA/oPool_validation_heatmap_ELISA.png): Figure 3B

### K<sub>D</sub> measurement via BLI
1. Compile BLI Kd raw data for plotting
``python3 script/validation/BLI/compile_BLI_result_batch.py``
    - Input files:
      - [./oPool_result/experimental_data/validation/BLI/Kd_measurements](./oPool_result/experimental_data/validation/BLI/Kd_measurements):raw data
    - Output files:
      - [./oPool_result/experimental_data/validation/BLI/Kd_compile](./oPool_result/experimental_data/validation/BLI/Kd_compile): compiled data

2. Plot kinetics data
``Rscript script/validation/BLI/plot_BLI_binding_batch.R``   
    - Input files:
      - [./oPool_result/experimental_data/validation/BLI/Kd_compile](./oPool_result/experimental_data/validation/BLI/Kd_compile): compiled data
    - Output files:
      - [./graph/validation/BLI/Kd_sensorgram](./graph/validation/BLI/Kd_sensorgram): all sensorgrams for kinetic measurements

3. Plot Kd heatmap
``Rscript script/validation/BLI/plot_Kd_heatmap.R``   
    - Output files:
      - [./graph/validation/BLI/H1_validation_heatmap.png](./graph/validation/BLI/H1_validation_heatmap.png) Figure 3D
      - [./graph/validation/BLI/H3_validation_heatmap.png](./graph/validation/BLI/H3_validation_heatmap.png) Figure 3D

### Competition validation via BLI
1. Plot competition results
``Rscript script/validation/BLI/plot_competition_validation.R``
    - Input files:
      - [./experimental_data/validation/BLI/oPool_competition_validation/oPool_competition_validation_sample_names.xlsx](./experimental_data/validation/BLI/oPool_competition_validation/oPool_competition_validation_sample_names.xlsx): Sample info for each file
      - [./experimental_data/validation/BLI/oPool_competition_validation/](./experimental_data/validation/BLI/oPool_competition_validation/) raw data
    - Output files:
      - [./graph/validation/BLI/oPool_competition_validation](./graph/validation/BLI/oPool_competition_validation): all sensorgrams for competition validation
      - [./experimental_data/validation/BLI/validated_antibody_competition_percentage.tsv](./experimental_data/validation/BLI/validated_antibody_competition_percentage.tsv)

2. Plot competition correlation
``Rscript script/validation/BLI/plot_competition_correlation.R``
    - Input files:
      - [./oPool_result/enrichment/competition_index.tsv](./oPool_result/enrichment/competition_index.tsv)    
      - [./experimental_data/validation/BLI/validated_antibody_competition_percentage.tsv](./experimental_data/validation/BLI/validated_antibody_competition_percentage.tsv)
    - Output files:
      - [./graph/validation/BLI/competition_index_vs_percentage.png](./graph/validation/BLI/competition_index_vs_percentage.png): Figure 4B

## Characterization of AG11-2F01 and 16.ND.92

### Structural analyses 
1. Plot structure overviews of AG11-2F01 and 16.ND.92
``pymol script/structural_analysis/overview.pml``   
    - Input files:
      - [./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb](./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb)
      - [./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb](./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/structural_analysis/PDB/2F01_overview.png](./graph/structural_analysis/PDB/2F01_overview.png): Figure 5A 
      - [./graph/structural_analysis/PDB/16ND_overview.png](./graph/structural_analysis/PDB/16ND_overview.png): Figure 5A 

2. Plot epitopes of AG11-2F01 and 16.ND.92
``pymol script/structural_analysis/epitope.pml``   
    - Input files:
      - [./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb](./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb)
      - [./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb](./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/structural_analysis/PDB/2F01_epitope.png](./graph/structural_analysis/PDB/2F01_epitope.png): Figure 5B
      - [./graph/structural_analysis/PDB/16ND_epitope.png](./graph/structural_analysis/PDB/16ND_epitope.png): Figure 5B

3. Plot CDRH3 interactions of AG11-2F01 and 16.ND.92
``pymol script/structural_analysis/interact_CDRH3.pml``   
    - Input files:
      - [./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb](./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb)
      - [./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb](./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/structural_analysis/PDB/2F01_CDRH3.png](./graph/structural_analysis/PDB/2F01_CDRH3.png): Figure 5C
      - [./graph/structural_analysis/PDB/16ND_CDRH3.png](./graph/structural_analysis/PDB/16ND_CDRH3.png): Figure 5C

4. Plot light chain interactions  of AG11-2F01 and 16.ND.92
``pymol script/structural_analysis/interact_LC.pml``   
    - Input files:
      - [./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb](./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb)
      - [./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb](./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb)
    - Output files:
      - [./graph/structural_analysis/PDB/2F01_LC.png](./graph/structural_analysis/PDB/2F01_LC.png): Figure 5D
      - [./graph/structural_analysis/PDB/16ND_LC.png](./graph/structural_analysis/PDB/16ND_LC.png): Figure 5D

5. Plot BSA of IGHD3-3 HA stem antibodies
``Rscript script/structural_analysis/plot_BSA_bar_chart.R``   
    - Input files:
      - [./experimental_data/structural_analysis/bsa_percentage.tsv](./experimental_data/structural_analysis/bsa_percentage.tsv)
    - Output files:
      - [./graph/structural_analysis/BSA_stacked_bar_chart.png](./graph/structural_analysis/BSA_stacked_bar_chart.png): Figure 5E

6. Plot IGHD3-3 contribution to VH paratopes
``Rscript script/structural_analysis/plot_3-3_BSA.R``   
    - Input files:
      - [./experimental_data/structural_analysis/3-3_BSA_percentage.tsv](./experimental_data/structural_analysis/3-3_BSA_percentage.tsv)
    - Output files:
      - [./graph/structural_analysis/BSA_IGHD3-3_bar_plot.png](./graph/structural_analysis/BSA_IGHD3-3_bar_plot.png): Figure 5F

7. Plot CDRH3 overlays of IGHD3-3 HA stem antibodies
``pymol script/structural_analysis/interact_LC.pml``   
    - Input files:
      - [./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb](./experimental_data/structural_analysis/PDB/SI06HA_2F01.pdb)
      - [./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb](./experimental_data/structural_analysis/PDB/SI06HA_16ND92.pdb)
      - [./experimental_data/structural_analysis/PDB/1G05_mono.pdb](./experimental_data/structural_analysis/PDB/1G05_mono.pdb)
      - [./experimental_data/structural_analysis/PDB/39.29_mono.pdb](./experimental_data/structural_analysis/PDB/39.29_mono.pdb)
      - [./experimental_data/structural_analysis/PDB/56.a.09_mono.pdb](./experimental_data/structural_analysis/PDB/56.a.09_mono.pdb)
      - [./experimental_data/structural_analysis/PDB/429_B01_mono.pdb](./experimental_data/structural_analysis/PDB/429_B01_mono.pdb)
      - [./experimental_data/structural_analysis/PDB/MEDI8852_mono.pdb](./experimental_data/structural_analysis/PDB/MEDI8852_mono.pdb)
      - [./experimental_data/structural_analysis/PDB/SIA28_mono.pdb](./experimental_data/structural_analysis/PDB/SIA28_mono.pdb)
    - Output files:
      - [./graph/structural_analysis/PDB/CDRH3_overlay.png](./graph/structural_analysis/PDB/CDRH3_overlay.png): Figure 5G

### Functional characterizations
8. Plot ELISA result heatmap
``Rscript script/functional_characterization/plot_EC50.R``   
    - Output files:
      - [./graph/functional_characterization/elisa_ec50_heatmap.png](./graph/functional_characterization/elisa_ec50_heatmap.png): Figure 6A

9. Plot micro-neutralization result heatmap
``Rscript script/functional_characterization/plot_IC50.R``   
    - Output files:
      - [./graph/functional_characterization/IC50_heatmap.png](./graph/functional_characterization/IC50_heatmap.png): Figure 6B

10. Plot in vivo experiment data
``Rscript script/functional_characterization/plot_in_vivo.R`` 
    - Input files:
      - [./experimental_data/functional_characterization/invivo_weight_loss.tsv](./experimental_data/functional_characterization/invivo_weight_loss.tsv)
      - [./experimental_data/functional_characterization/invivo_survival.tsv](./experimental_data/functional_characterization/invivo_survival.tsv)
      - [./experimental_data/functional_characterization/invivo_lung_titer.tsv](./experimental_data/functional_characterization/invivo_lung_titer.tsv)
    - Output files:
      - [./graph/functional_characterization/invivo_weight_loss_2F01.png](./graph/functional_characterization/invivo_weight_loss_2F01.png): Figure 6C
      - [./graph/functional_characterization/invivo_weight_loss_16ND92.png](./graph/functional_characterization/invivo_weight_loss_16ND92.png): Figure 6D
      - [./graph/functional_characterization/invivo_survival_2F01.png](./graph/functional_characterization/invivo_survival_2F01.png): Figure 6E
      - [./graph/functional_characterization/invivo_survival_16ND92.png](./graph/functional_characterization/invivo_survival_16ND92.png): Figure 6F
      - [./graph/functional_characterization/invitro_lung_titer_Vero_2F01.png](./graph/functional_characterization/invitro_lung_titer_Vero_2F01.png): Figure 6G
      - [./graph/functional_characterization/invitro_lung_titer_Vero_16ND92.png](./graph/functional_characterization/invitro_lung_titer_Vero_16ND92.png): Figure 6H


