#!/usr/bin/env bash
# This script documents the reproducible metagenome bioinformatics pipeline by Jiang Lab
# Author: Sang Yimeng
# Time: 2023/03/30

#Need software list (recommended to install everything using conda)
conda install -c bioconda seqkit
conda install -c bioconda emboss
conda install -c bioconda fastqc
conda install -c bioconda trimmomatic
conda install -c bioconda megahit
conda install -c bioconda prodigal
conda install -c bioconda cd-hit
conda install -c bioconda salmon
#Refer to github (https://github.com/bxlab/metaWRAP) for Metawrap installation.

# set parameter
threads=60 #set maximum number of threads
workspace=~/workspace/
1rawdatap=~/workspace/rawdata/
2qcresult=~/workspace/qcresult/
3cleandatap=~/workspace/cleandata/
4assemblyp=~/workspace/assembly/
5genep=~/workspace/gene/
6abundancep=~/workspace/abundance
7binningp=~/workspace/binning
mkdir ${1rawdatap} ${2qcresult} ${3cleandatap} ${4assemblyp} ${5genep} ${6abundancep} ${7binningp}

#Step1: Use Fastqc assess data quality and use Trimmomatic filterted
# Check the data quality
cd ${workspace}
for i in `ls ${1rawdatap}`
do
	echo "fastqc -o ${2qcresult} ${i}" >> do.sh
done
nohup parallel -j ${threads} < do.sh &
rm do.sh
# Cut out low quality sequences
for i in `ls ${1rawdatap}/*_1.fq.gz`
do
	filename=echo ${i}|awk -F '[/]' '{print $NF}'|cut -d "_" -f1
	trimmomatic PE -threads ${threads} -phred33 ${1rawdatap}/${filename}_1.fq.gz ${1rawdatap}/${filename}_2.fq.gz ${3cleandatap}/${filename}_1_pairer.fq.gz ${3cleandatap}/${filename}_1_unpairer.fq.gz ${3cleandatap}/${filename}_2_pairer.fq.gz ${3cleandatap}/${filename}_2_unpairer.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
# filename_1_pairer.fq.gz is both terminal sequences are preserved
# filename_1_unpairer.fq.gz is keeping only one end of the sequence
# Recheck data qyality
for i in `ls ${3cleandatap}`
do
        echo "fastqc -o ${2qcresult} ${i}" >> do.sh
done
nohup parallel -j ${threads} < do.sh &
rm do.sh

#Step2: assembly, Gene prediction, cluster, quantitfy
# assembly megahit (Attention: need long time and high server load)
cd ${3cleandatap}
for i in `ls *_1_pairer.fq.gz`
do
	filename=echo ${i}|cut -d "_" -f1
	megahit -1 ${filename}_1_pairer.fq.gz -2 ${filename}_2_pairer.fq.gz -o ${4assemblyp}/${filename}/
	seqkit stat ${4assemblyp}/${filename}/final.contigs.fa #check the result status
done
# Prediction metaProdigal
cd ${4assemblyp}
for i in `ls ${3cleandatap}/*_1_pairer.fq.gz`
do
	filename=echo ${i}|awk -F '[/]' '{print $NF}'|cut -d "_" -f1
	prodigal -i ${filename}/final.contigs.fa  -d ${5genep}/${filename}.fa -o ${5genep}/${filename}.gff -p meta -f gff > ${5genep}/prodigal.log 2>&1
	grep -c '>' ${5genep}/${filename}.fa #Statistical gene number
	grep -c 'partial=00' ${5genep}/${filename}.fa #Count the number of complete genes
	grep 'partial=00' ${5genep}/${filename}.fa | cut -f1 -d ' '| sed 's/>//' > ${5genep}/${filename}_full_length.id
	seqkit grep -f ${5genep}/${filename}_full_length.id ${5genep}/${filename}.fa > ${5genep}/${filename}_full_length.fa
	seqkit stat ${5genep}/${filename}_full_length.fa #Extraction of complete genes (such as rings of bacterial genomes)
done
# Clustercd (Attention: need long time and high server load)
cd ${5genep}
for i in `ls ${3cleandatap}/*_1_pairer.fq.gz`
do
	filename=echo ${i}|awk -F '[/]' '{print $NF}'|cut -d "_" -f1
	time cd-hit-est -i ${filename}.fa -o ${filename}_nucleotide.fa -aS 0.9 -c 0.95 -G 0 -g 0 -T ${threads} -M 0
done
cat *_nucleotide.fa >> allgene.fa
time cd-hit-est -i ${5genep}/allgene.fa -o ${5genep}/unigene.fa -aS 0.9 -c 0.95 -G 0 -g 0 -T ${threads} -M 0
transeq -sequence ${5genep}/unigene.fa -outseq ${5genep}/uniprotein.fa -trim Y # translation nucleotide to protein
# quantitfy salmon
cd ${6abundancep}
time salmon index -t ${5genep}/unigene.fa -p ${threads} -i ${6abundancep}/index
for i in `ls ${3cleandatap}/*_1_pairer.fq.gz`
do
	filename=echo ${i}|awk -F '[/]' '{print $NF}'|cut -d "_" -f1
	echo "salmon quant -i ${6abundancep}/index -1 ${3cleandatap}/${filename}_1_pairer.fa -2 ${3cleandatap}/${filename}_2_pairer.fa -o ${6abundancep}/${filename}.quant" >> do.sh
done
nohup parallel -j ${threads} < do.sh &
rm do.sh
salmon quantmerge --quants ${6abundancep}/*.quant -o ${6abundancep}/gene.TPM
salmon quantmerge --quants ${6abundancep}/*.quant --column NumReads -o ${6abundancep}/gene.count
# TPM adjust result in gene.TPM as relative abundance

#Step3: Binning
# Using MetaWRAP(https://github.com/bxlab/metaWRAP)
# Attention: need long time and high server load
# Usage tutorial in: https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md
cd ${7binningp}
nohup metawrap binning -o ${7binningp} -t ${threads} -a ${5genep}/*_full_length.fa --metabat2 --maxbin2 --concoct ${3cleandatap}/*pairer.fq.gz &

# After this pipeline, we get reads, assemblies and binnings for our sample.
# These three approaches (read profiling, assembly profiling, and binning profiling) can be used together or separately to gain a comprehensive understanding of the taxonomic and functional composition of a microbial community in a metagenomic sample.
# 1.Read profiling can be performed using tools such as Kraken, Kaiju, and MetaPhlAn
# 2.Assembly profiling and binning profiling might annotation in different way such as blast, hmmsearch and so on.
# Downstream analysis requires thoughtful experimental design and strategy and cannot be shown by a single script.
