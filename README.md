# oPool<sup>+</sup> display: A high-throughput cell-free screening platform for natively paired antibodies
Welcome! This README describes the workflow in the manuscript: High-throughput synthesis and specificity characterization of natively paired antibodies using oPool<sup>+</sup> display (Ouyang et al., 2025)

## Contents
- [Introduction](#introduction)
- [Environment setup](#environment-setup)
- [oPool<sup>+</sup> display design](#oPool<sup>+</sup>-display-design)
- [oPool<sup>+</sup> display results](#oPool<sup>+</sup>-display-results)


## Introduction
oPool<sup>+</sup> display combines oligo pool synthesis and mRNA display to construct and characterize the specificity of many natively paired antibodies in parallel. As a proof-of-concept, we applied oPool<sup>+</sup> display to rapidly screen the binding activity of >300 previously uncharacterized influenza hemagglutinin (HA) antibodies against 9 HA variants via 16 different screens. Over 5,000 binding tests were performed in 3-5 days. This repository, therefore, contains two parts: 

    1. the [oPool_design](oPool_design/) folder: oligo sequence designs for library assembly; 
    2. the [oPool_display](oPool_display/) folder: screening/validation results and analyses.

## Environment setup 
To create the same working environment:
```bash
conda env create -f environment.yml
```
## Oligo pool seqeunce design

Please proceed to the [oPool_design](oPool_design/) folder

## Screening and experimental data analyses

Please proceed to the [oPool_display](oPool_display/) folder