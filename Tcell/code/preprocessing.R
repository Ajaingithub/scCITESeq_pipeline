#### This script is scCITEseq Ines Run01 Demultiplexing using Cell Ranger regressing unwanted Effect and integrating Modalities and samples keeping Feature 4000. Also instead of scTransform using Normalize Data #

#### Before Command ########
# https://satijalab.org/seurat/articles/visualization_vignette.html
# Download the fastq data from the Novogene on the screen. This will take 2-2.5 days. Since data is not demultiplexed. Try to download it in parallel but keep track everything is downloaded
# Check the bytes as well as MD5 for your files using md5sum filename with the MD5.txt they have provided
# ls *.fastq.gz | awk '{print "md5sum "$0" >> checking_MD5.txt"}' > MD5_check.sh
# Merge the data for both RNA and ADT. all R1, R2, and I1. Since we are using Antibody tagged oligos, GEX. We will get
# 2 fastqs. 1. GEX, 2. CSP. Need to make a new folder with the fastq from CSP. Remove the first read from all R1, R2, and I1
# Read this https://kb.10xgenomics.com/hc/en-us/articles/4407386498957-I-used-antibody-tags-for-cell-surface-protein-capture-and-cell-hashing-with-Single-Cell-3-chemistry-How-can-I-use-Cell-Ranger-to-analyze-my-data-
# removing 1st four lines from the CSP from R1, R2, and I1. In order to have separate identity from the CSP for multiplex capture.
# zcat Exp01_CSP_combined_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3_S1_L001_I1_001.fastq.gz | sed '1,4d' > Exp01_CSP_combined_del_4_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3_S1_L001_I1_001.fastq
# ls *.fastq.gz | sed 's/_CSP_/\t/g' | awk '{print "zcat "$1"_CSP_"$2" | sed \x271,4d\x27 > /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run02_and_Run03/combined_raw_data/Exp02_and_03_CSP_multiplex_capture/"$1"_CSP_4_del"$2" &"}' > multiplex_fastq.sh
# please make change \x271 has put q make it '1 in all the files
# Also remove the gz at the end of the file
# gzip Exp01_CSP_combined_del_4_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3_S1_L001_I1_001.fastq
# ls *.fastq | awk '{print "gzip "$0" &"}' > zipping.sh
# Now you have three folder 1. GEX, 2. CSP, 3. multiplexing capture

# Need to make CMO reference for demultiplexing your sample
# id,name,read,pattern,sequence,feature_type
# HTO1,HTO1,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GTCAACTCTTTAGCG,Multiplexing Capture
# HTO2,HTO2,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGATGGCCTATTGGG,Multiplexing Capture
# HTO3,HTO3,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TTCCGCCTCTCTTTG,Multiplexing Capture
# HTO4,HTO4,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AGTAAGTTCAGCGTA,Multiplexing Capture
# HTO5,HTO5,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AAGTATCGTTTCGCA,Multiplexing Capture
# HTO6,HTO6,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GGTTGCCAGATGTCA,Multiplexing Capture
# HTO7,HTO7,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGTCTTTCCTGCCAG,Multiplexing Capture
# HTO8,HTO8,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CTCCTCTGCAATTAC,Multiplexing Capture
# HTO9,HTO9,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CAGTAGTCACGGTCA,Multiplexing Capture
# HTO10,HTO10,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,ATTGACCCGCGTTAG,Multiplexing Capture

# Also need to make a protein reference keep the negative control also to find if any cluster has high negative control in that
# id,name,read,pattern,sequence,feature_type
# CD4_protein_Total_SeqB,CD4_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGTTCCCGCTCAACT,Antibody Capture
# CD8a_protein_Total_SeqB,CD8a_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GCTGCGCTTTCCATT,Antibody Capture
# CD45RA_protein_Total_SeqB,CD45RA_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TCAATCCTTCCGCTT,Antibody Capture
# CD45RO_protein_Total_SeqB,CD45RO_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CTCCGAATCATGTTG,Antibody Capture
# CD197_CCR7_protein_Total_SeqB,CD197_CCR7_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AGTTCAGTCAACCGA,Antibody Capture
# CD95_Fas_protein_Total_SeqB,CD95_Fas_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CCAGCTCATTAGAGC,Antibody Capture

# Make the multiconfig.csv file, This include 
# [gene-expression]
# reference,/diazlab/data3/.abhinav/.immune/resources/10x_reference/refdata-gex-GRCh38-2024-A/
# cmo-set,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/analysis/preprocessing/required_files/Run01_CMO_reference.csv
# no-secondary,TRUE
# no-bam,TRUE

# [feature]
# reference,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/analysis/preprocessing/required_files/Total_seqB_with_HTO_reference.csv

# [libraries]
# fastq_id,fastqs,lanes,feature_types,subsample_rate
# 428117-GEX_1-Z0193-CTGCACATTGTAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/GEX/,any,Gene Expression,
# 428720-GEX_1-Z0193-CTGCACATTGTAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/GEX/,any,Gene Expression,
# 428744-GEX_1-Z0193-CTGCACATTGTAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/GEX/,any,Gene Expression,
# 428117-CITE_1-Z0097-CAGTCAGTTGCAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/CITE/,any,Antibody Capture,
# 428720-CITE_1-Z0097-CAGTCAGTTGCAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/CITE/,any,Antibody Capture,
# 428744-CITE_1-Z0097-CAGTCAGTTGCAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/CITE/,any,Antibody Capture,
# 428117-CITE_4_del1-Z0097-CAGTCAGTTGCAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/HTO/,any,Multiplexing Capture
# 428720-CITE_4_del1-Z0097-CAGTCAGTTGCAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/HTO/,any,Multiplexing Capture
# 428744-CITE_4_del1-Z0097-CAGTCAGTTGCAGAT,/diazlab/data3/.abhinav/.immune/Mayo/Tcells/rawdata/HTO/,any,Multiplexing Capture

# [samples]
# sample_id,cmo_ids,description
# T_Y1_33M_Proliferating-Control,HTO1,T_Y1_33M_Proliferating-Control
# T_Y2_30M_Proliferating-Control,HTO2,T_Y2_30M_Proliferating-Control
# T_Y3_30M_Proliferating-Control,HTO3,T_Y3_30M_Proliferating-Control
# T_Y4_22F_Proliferating-Control,HTO4,T_Y4_22F_Proliferating-Control
# T_Y5_30F_Proliferating-Control,HTO5,T_Y5_30F_Proliferating-Control
# T_O1_71M_Proliferating-Control,HTO6,T_O1_71M_Proliferating-Control
# T_O2_79M_Proliferating-Control,HTO7,T_O2_79M_Proliferating-Control
# T_O3_77M_Proliferating-Control,HTO8,T_O3_77M_Proliferating-Control
# T_O4_72F_Proliferating-Control,HTO9,T_O4_72F_Proliferating-Control
# T_O5_74F_Proliferating-Control,HTO10,T_O5_74F_Proliferating-Control

# source ~/.bashrc
# module load cellranger/9.0.1
# cellranger multi --id=Run01_CR9 --csv=/diazlab/data3/.abhinav/.immune/Mayo/Tcells/analysis/preprocessing/required_files/Run01_multiconfig_CR9.csv

### Running it on Slurm
# #!/bin/bash
# #SBATCH --job-name=Run02
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=8
# #SBATCH --mem=200G
# #SBATCH --time=4-00:00:00
# #SBATCH --output=./log/Run02_CR9_%j.out
# #SBATCH --error=./log/Run02_CR9_%j.err
# #SBATCH --mail-user=abhinav.jain@ucsf.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50

# source ~/.bashrc
# ml load cellranger/9.0.1

# cellranger multi \
#     --id=Run02_CR9 \
#     --csv=/diazlab/data3/.abhinav/.immune/Mayo/Tcells/analysis/preprocessing/required_files/Run02_multiconfig_CR9.csv \
#     --localcores=8 \
#     --localmem=200
