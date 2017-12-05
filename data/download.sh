#!/bin/bash

#=======================================================================
# this script is suppost do document all downlaods in the data folder
#=======================================================================

# set some variables here:
BIN=../bin
mkdir -p ${BIN}

#=======================================================================
# Hi-C data from Rao et al 2014 Cell
#=======================================================================
mkdir -p Rao2014

# K562 inter-chromosomal matrices (~6.4 GB in .gz)
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FK562%5Finterchromosomal%5Fcontact%5Fmatrices%2Etar%2Egz

# K562 inta-chromosomal matrices (~6.4 GB in .gz)
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FK562%5Fintrachromosomal%5Fcontact%5Fmatrices%2Etar%2Egz

# GM12878 intra and inter chromosomal contact data:
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Finterchromosomal%5Fcontact%5Fmatrices%2Etar%2Egz
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Fintrachromosomal%5Fcontact%5Fmatrices%2Etar%2Egz

#~ tar xvfz Rao2014/*.tar.gz

# unzip all
#gunzip Rao2014/*.gz
tar xvfz Rao2014/GSE63525_K562_interchromosomal_contact_matrices.tar.gz -C Rao2014
tar xvfz Rao2014/GSE63525_K562_intrachromosomal_contact_matrices.tar.gz -C Rao2014
tar xvfz Rao2014/GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz -C Rao2014
tar xvfz Rao2014/GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz -C Rao2014

RAO_CELLS="GM12878_primary+replicate HMEC HUVEC HeLa IMR90 K562 KBM7 NHEK"

for CELL in ${RAO_CELLS} ; do
    # download
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_HiCCUPS_looplist.txt.gz
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz

    # unzip 
    gunzip Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
    gunzip Rao2014/GSE63525_${CELL}_HiCCUPS_looplist.txt.gz
    gunzip Rao2014/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz
        
    # re-format TADs into bed file
    tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt |cut -f 1-3 | sed -e 's/^/chr/' > Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 
done

#=======================================================================
# Capture Hi-C data from Mifsud2015
#=======================================================================
mkdir -p Mifsud2015
wget -P Mifsud2015 http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip
unzip Mifsud2015/E-MTAB-2323.additional.1.zip -d Mifsud2015
