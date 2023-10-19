#!/usr/local/bin/Rscript
# R 4.1.0-icelake HPC
#---------------------------------------------------------------------------------
# SLX-21013 (GTC251), GADD34 KO/WT with different time post irradiation
# non-irradiated, 2wk, 4wk and 7wk; Mouse 
# each type has 3 replicates.
# CAD_jcg23_0001; Ensemble GRCm39
# Information about your samples related to bioinformatics:
# 1. These were Illumina TruSeq Stranded mRNA libraries.
# 2. PE50 was sequenced on Illumina NovaSeq6000 platform.
# 3. The average library size (including adapters) was 340 bp.
# 4. Index type was Truseq RNA CD.
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CAD-BFX/CAD_jcg23_0001
#
#
# Analysis Performed by Xiaohui Zhao
# Department of medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#---------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+---- Basic settings and libraries calling  ---------+")

## module load R/4.1.0-icelake

setwd("/home/xz289/rds/rds-cad_bfx_projects-OsGdhNr6c2U/CAD-Projects/jcg23/CAD_jcg23_0001/Original_Data/SLX-21013")
library("readxl")

cruk <- read.csv("SLX-21013.HMHHWDRXY.s_1.contents.csv", header=T)
samT <- read_excel("SLX-21013_Summary_Table.xlsx")
colnames(samT) <- c("Sample.name", "Condition", "WeeksPI")
samTable <- merge(samT, cruk, by = "Sample.name")
samTable <- samTable[order(samTable$Barcode), ]
samTable$time <- rep(c("0w", "2w", "4w", "7w"), each=6)
samTable$rep  <- c(rep(c(1:3), length=6), c(1:4), c(1:2), rep(c(1:3), length=12))
fastq1   <- list.files(path = ".", pattern = "*.r_1.fq.gz")
fastq2   <- list.files(path = ".", pattern = "*.r_2.fq.gz")
samTable$fastq_1 <- fastq1
samTable$fastq_2 <- fastq2
samTable$sample  <- paste0(samTable$Condition, "_T", samTable$time, "_REP", samTable$rep)
samTable$strandedness <- "reverse"

write.csv(samTable, file = "CAD_jcg23_0001-SampleTable.csv", row.names = F)
write.csv(samTable[,c("sample", "fastq_1", "fastq_2", "strandedness")],
          file = "CAD_jcg23_0001-nextflow_SampleTable.csv", row.names=F)

##---------------- FINISH----------------------------------------##


