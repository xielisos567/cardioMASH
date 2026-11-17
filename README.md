## Introduction
cardioMASH is designed to determine the cardiovascular disease associated metabolic gene clusters in a set of given genomic data. The input files are expected to involve the bacterial genomes, the species-level genome bins (SGBs), the PacBio long reads or metagenome-assembled genomes (MAGs).

**Note: The core functionality of cardioMASH is primarily designed for BLASTN (nucleotide-to-nucleotide alignment), which is highly recommended for analyzing bacterial genomes and SGBs. The TBLASTN and BLASTX modes are provided as auxiliary functions for metagenomic data or annotated CDS.

## System requirements
### Operating System
* Linux (e.g., Ubuntu 20.04 or CentOS 7) or macOS.
* The software has been tested on Ubuntu 20.04 LTS.

### Dependencies
* [BLAST: tested with v2.16.0+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [Prodigal: tested with v2.6.3](https://github.com/hyattpd/Prodigal)
* [BBMap: tested with version 39.26](https://github.com/BioInfoTools/BBMap)
* [DIAMOND: tested with v0.8.36.98](https://github.com/bbuchfink/diamond/releases)
* [Python>=3.7,<3.12](https://www.python.org/downloads/)
* Python modules: 
    * [Pandas < 1.4.0](http://pandas.pydata.org/pandas-docs/stable/install.html)
    * [numpy](https://numpy.org/)

### Hardware requirements
 *Standard hardware: No non-standard hardware is required. Runs on standard Linux/Unix systems.
 
 *Recommended configuration: CPU threads ≥ 8, RAM ≥ 64 GB.

## Installation guide

#### Step 1. Download databases
Sequence alignment requires the [cvdMGC database](https://github.com/xielisos567/cvdMGC) before usage of the cardioMASH. You can download the database (Redundant_CVD_associated_genes) manually. 

#### Step 2. Download main scripts
git clone https://github.com/xielisos567/cardioMASH

cd cardioMASH/core_functionality

#### Step 3. Install the Dependencies
We recommend using [Mamba](https://github.com/mamba-org/mamba) to install all dependencies automatically. 
An `environment.yml` file is provided in the repository. You can create the environment by running:
```
mamba env create -f environment.yml
mamba activate cardioMASH
```
Typical install time: Approximately 5-10 minutes on a normal desktop computer (excluding database download time).

## Demo
The demo dataset (15_genomes/), the genome list (15_genomes.txt), the species-specific average nucleotide identity file (ANI_result.txt), and the target gene label file (gene_lable) are all provided in the repository.
 * Expected run time: With the core functionality (BLASTN), less than 5 minutes on a normal desktop computer (8 threads, 16GB RAM).
 * Expected output: The demo will generate a combined alignment result (combine_blastn_results_all.out), a filtered result based on species-specific ANI (15_genomes_ssANI.txt), and the final annotation output in the final_result2/ directory, which summarizes the detected MGCs and associated taxonomy.

### Instructions to run on demo data

#### Step 1: Batch build BLAST database for 15 genomes
```
for i in `ls 15_genomes/`
do
  makeblastdb -in ./15_genomes/$i -dbtype nucl
done
```
**Note: Before sequence alignment with cvdMGC, you need to download Redundant_CVD_associated_genes.zip from cvdMGC repository, and merge all CDS into a single combine_CVD_cds.fasta file. Additionally, use the sed command to replace spaces with underscores:
```
sed -i 's/ /_/g' combine_CVD_cds.fasta
```
#### Step 2: Batch sequence alignment with the cvdMGC database
```
for i in `cat ./15_genomes.txt`
do
  blastn -query ../cvdMGC/combine_CVD_cds.fasta -db ./15_genomes/${i}.fna -evalue 1e-5 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 1 -out result/${i}.bls.fmt6 -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 6
done
```
#### Step 3: Merge sequence alignment results
```
cd result
for i in `ls ./`
do 
  cat ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../combine_blastn_results_all.out
done
cd ..
```
#### Step 4: ssANI filtering
```
python3 data-filter_ssANI2.py ANI_result.txt combine_blastn_results_all.out 15_genomes_ssANI.txt
```
#### Step 5: cardioMASH annotation
(The gene_lable file contains the target gene clusters and associated genes to be extracted)
```
python3 gene_cluster_extract_final_ssANI6-fragment.py gene_lable 15_genomes_ssANI.txt final_result2/15_genomes_annotation
```
**Expected output structure of `final_result2/`:
The `final_result2/` directory contains the final annotation outputs. The directory structure is as follows:

```
final_result2/
├── 15_genomes_annotation.complete.rmdup.xls
├── 15_genomes_annotation.complete.stat.xlsx
├── 15_genomes_annotation.fail.xls
├── 15_genomes_annotation.filter.xls
├── 15_genomes_annotation.fragmentary.rmdup.xls
├── 15_genomes_annotation.fragmentary.stat.xlsx
├── 15_genomes_annotation.fragmentary.xls
├── 15_genomes_annotation.gene.data.xls
└── 15_genomes_annotation.gene.stat.xlsx
```

## Instructions for use

### Core Function: Running BLASTN
(Recommended for bacterial genomes or SGBs)

#### Step 1. BLASTN through the nucleotide database of cvdMGC 
(Note: Before running BLASTN, concatenate the gene IDs in the headers of the merged CDS sequences using sed, replacing spaces with underscores, e.g.: sed -i 's/ /_/g' all_redundant_genes.fasta)
```
blastn -query cvdMGC/all_redundant_genes.fasta -db <makeblastdb_genomes> -evalue 1e-5 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 1 -out result/<BLASTN_alignment_results> -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 6
```
#### Step 2: Merge sequence alignment results
```
cd result
for i in `ls ./`
do 
  cat ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../<BLASTN_alignment_results_merged>
done
cd ..
```
#### Step 3. Filtering the aligned results through ssANI threshold
```
python3 data-filter_ssANI2.py ANI_result.txt <BLASTN_alignment_results_merged> <output_file_ssANI>
```
#### Step 4. Architecture search physically clustered genes and sequence extraction
```
python3 gene_cluster_extract_final_ssANI6-fragment.py gene_lable <output_file_ssANI> <output_result_ssANI_statics>
```

### Auxiliary Function1: Running BLASTX
(For metagenome assembly and annotated CDS)
#### Step 0. makedb:
```
cd cvdMGC/
diamond makedb --in all_redundant_proteins.faa --db all_redundant_proteins --threads 8
```
#### Step 1. BLASTX through the protein database
```
diamond blastx --db cvdMGC/all_redundant_proteins --query <annotated_cds_seqs> --out result/<BLASTX_alignment_results> --more-sensitive --evalue 1e-5 --id 40 --outfmt 6 --query-cover 80
```
#### Step 2. Merge sequence alignment results
```
cd result
for i in `ls ./`
do 
  cat ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../<BLASTX_alignment_results>
done
```
#### Step 3. Data format conversion
```
python3 data_format_trans.py <BLASTX_alignment_results> <output_trans_blastx> c
```
#### Step 4. Calculate the final results
```
python3 gene_cluster_extract_blastx.py gene_lable <output_trans_blastx> final_result2/<output_result_blastx_statics>
```

### Auxiliary Function2: Running TBLASTN
(For protein-to-genome alignment)
#### Step 0. makedb:
```
cd genomes/
for i in `ls ./`
do
  makeblastdb -in ./$i -dbtype nucl
done
```
#### Step 1. TBLASTN through the protein database of cvdMGC
```
tblastn -query cvdMGC/all_redundant_proteins.faa -db <makeblastdb_genomes> -evalue 1e-5 -qcov_hsp_perc 80 -max_target_seqs 1 -out result/<TBLASTN_alignment_results> -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 16
```
#### Step 2. Merge sequence alignment results
```
cd result
for i in `ls ./`
do 
  cat ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../<TBLASTN_alignment_results>
done
```
#### Step 3. Data format conversion
```
python3 data_format_trans.py <TBLASTN_alignment_results> <output_trans_tblastn> s
```
#### Step 4. Calculate the final results
```
python3 gene_cluster_extract_tblastn.py gene_lable <output_trans_tblastn> final_result2/<output_result_tblastn_statics> tblastn
```


## Copyright
Shulei Jia, E-mail: jiashu320lei@126.com; jiashulei@tmu.edu.cn

School of Basic Medical Sciences, Tianjin Medical University, Tianjin, 300070, China

Institute of Microbiology, Chinese Academy of Sciences, Beijing, 100101, China