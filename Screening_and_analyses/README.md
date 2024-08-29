# Screening results and data analyses

## Input files
* [./ref_files/300lib_Abs.csv](./ref_files/300lib_Abs.csv): Table S1, information of selected antibodies
* [./ref_files/300lib.tsv](./ref_files/300lib.tsv): Reference sequences of the natively paired antibody design
* [./ref_files/neg_abs_list.tsv](./ref_files/neg_abs_list.tsv): List of the 25 HA head antibodies (negative controls)
* [./ref_files/sample_name.tsv](./ref_files/sample_name.tsv): Sample names of PacBio sequencing files
* Raw read (PacBio CCS) files in fastq format from NIH SRA database [BioProject PRJNA1150188](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1150188)
* [./data](./data/): Experimental data (validation, structural & functional characterization)

## Analyses of mRNA display results

### Data preparation
1. Filter the CCS reads based on quality and ROI, then trim the adaptors   
``python3 script/filter_from_fastq.py``
    - Input files: PacBio CCS reads
    - Output files will be placed in a folder named fastq_filtered/

### Counting reads, identify scFvs, and calculate enrichment for both replicates
2. Counting unique reads, rename samples
``python3 script/fastq2count.py``   
    - Input files:
      - Merged read files in fastq_merged/ folder
      - [./ref_files/sample_name.tsv](./ref_files/sample_name.tsv)
    - Output files:
      - [./result/PacBio/mut_nuc_count.tsv](./result/PacBio/mut_nuc_count.tsv)

3. Identify natively paired scFvs with no mutations, then calculate frequncy and enrichment
``python3 script/identify_scfv_300lib.py``   
    - Input files:
      - [./result/PacBio/mut_nuc_count.tsv](./result/PacBio/mut_nuc_count.tsv)
      - [./ref_files/300lib_Abs.csv](./ref_files/300lib_Abs.csv)
      - [./ref_files/300lib.tsv](./ref_files/300lib.tsv)
    - Output files:
      - [./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv](./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv)
      - [./result/PacBio/oPool_screen_enrichment.tsv](./result/PacBio/oPool_screen_enrichment.tsv)

### Assembly assessment
4. Plot scFv frequency correlation between two assembly replicates
``Rscript script/plot_QC_input.R``   
    - Input files:
      - [./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv](./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv)
    - Output files:
      - [./graph/input_QC.png](./graph/input_QC.png)

### Screening result analyses
5. Plot H1 stem screen correlation between replicates
``Rscript script/plot_QC_H1_screen.R``   
    - Input files:
      - [./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv](./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv)
    - Output files:
      - [./graph/H1_screen_correlation.png](./graph/H1_screen_correlation.png)

6. Plot H3 stem screen correlation between replicates
``Rscript script/plot_QC_H3_screen.R``   
    - Input files:
      - [./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv](./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv)
    - Output files:
      - [./graph/H3_screen_correlation.png](./graph/H3_screen_correlation.png)

7. Plot H1 stem screen enrichment vs. H1 stem screen enrichment 
``Rscript script/plot_QC_H1_H3.R``   
    - Input files:
      - [./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv](./result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv)
    - Output files:
      - [./graph/H1_vs_H3.png](./graph/H1_vs_H3.png)

## Analyses of experimental validation results

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


