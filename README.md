## Introduction
cardioMASH is designed to determine the cardiovascular disease associated metabolic gene clusters in a set of given genomic data. The input files are expected to involve the bacterial genomes, the species-level genome bins(SGBs), the PacBio long reads or metagenome-assembled genomes (MAGs).

## Dependencies
* [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [BBMap](https://github.com/BioInfoTools/BBMap)
* [Python3](https://www.python.org/downloads/)
* Python modules: 
    * [Pandas < 1.4.0](http://pandas.pydata.org/pandas-docs/stable/install.html)
    * [numpy](https://numpy.org/)

## Installation
Recommended configuration:  
CPU threads ≥ 8  
RAM ≥ 64 Gb

### Repository
#### Step 1. Download databases
Sequence alignment requires [cvdMGC database](http://43.143.223.231:88/CVD_DATA_database/) before usage of the cardioMASH.
You can download the database manually. The sequence alignment was conducted between query and the cvdMGC database at designed threshold.

#### Step 2. Download main scripts
`git clone https://github.com/xielisos567/cardioMASH`

#### Step 3. Install the Dependencies
The [Dependencies](#Dependencies) are required to be installed and added to the system `$PATH`

## Running
### 1. Running BLASTN
#### Step 1. BLASTN through the nucleotide database of cvdMGC(alignment between nucleotide sequences and genome)
```
blastn -query cvdMGC/all_cds.fna -db <makeblastdb_genomes> -evalue 1e-5 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 1 -out <BLASTN_alignment_results> -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 4
```

#### Step 2. Filtering the aligned results through species-specific average nucleotide identity (ssANI) threshold
```
python3  data-filter_ssANI2.py  ANI_result.txt  <BLASTN_alignment_results>  <output_file_ssANI>
```

#### Step 3. Architecture search physically clustered genes and sequence extraction
```
python3  gene_cluster_extract_final_ssANI5-fragment  gene_lable  <output_file_ssANI>  <output_result_ssANI_statics>
```

##### Run examples
See examples in the folder of blastn
```
URL: https://github.com/xielisos567/cardioMASH/blastn
```

## 2. Running TBLASTN
#### Step 1. TBLASTN through the protein database of cvdMGC(alignment between protein and genome)
```
tblastn -query cvdMGC/all_translated_proteins.faa -db <makeblastdb_genomes> -evalue 1e-5 -qcov_hsp_perc 80 -max_target_seqs 1 -out <TBLASTN_alignment_results> -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 4
```

#### Step 2. Data format conversion (letter s represents the format of genome number in the first column and gene number in the second column)
```
python3 data_format_trans.py <TBLASTN_alignment_results> <output_trans_tblastn> s
```

#### Step 3. Calculate the final results
```
python3 gene_cluster_extract_tblastn.py gene_lable <output_trans_tblastn> <output_result_tblastn_statics> tblastn
```

#### Run examples
See examples in the folder of tblastn_blastx
```
URL: https://github.com/xielisos567/cardioMASH/tblastn_blastx/tblastn
```

## 3. Running BLASTX
#### Step 1. BLASTX through the protein database of cvdMGC(alignment between protein and cds annotated with prodigal):
```
time diamond blastx --db cvdMGC/all_translated_proteins.faa --query <annotated_cds_seqs> --out <BLASTX_alignment_results> --more-sensitive --evalue 1e-5 --id 40 --outfmt 6 --query-cover 80
```

#### Step 2. Data format conversion (letter c represents the format of the second column with both genome and gene numbers)
```
python3 data_format_trans.py <BLASTX_alignment_results> <output_trans_blastx> c
```

#### Step 3. Calculate the final results
```
python3 gene_cluster_extract_blastx.py gene_lable <output_trans_blastx> <output_result_blastx_statics>
```

#### Run examples
See examples in the folder of tblastn_blastx.
```
URL: https://github.com/xielisos567/cardioMASH/tblastn_blastx/blastx
```

## Recommend && Notice
```
Genomic analysis of genomes from bacterial strains or SGBs (metagenome binning), BLASTN is recommended for sequence alignments
```
```
Genomic analysis of metagenome such as the assembled genomes, BLASTX is recommended for sequence alignments
```


## Copyright
Shulei Jia, E-mail: jiashulei@tmu.edu.cn  
School of Basic Medical Sciences, Tianjin Medical University, Tianjin, 300070, China
Institute of Microbiology, Chinese Academy of Sciences, Beijing, 100101, China
