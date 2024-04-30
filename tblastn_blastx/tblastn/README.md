###tblastn: sequence alignment of the protein sequence to the genome database, the results are filtered through the pipeline of ssANI threshold:
```
for i in `ls 10_genomes_db/`
do
	tblastn -query combine_all_caiTABCDE_length_filter.faa -db 10_genomes_db/$i -evalue 1e-5 -qcov_hsp_perc 80 -max_target_seqs 1 -out blast_result_final/${i}_10.bls.fmt6 -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 4
done
```

###Step 1, data format conversion(letter s represents the format of genome number in the first column and gene number in the second column):
```
python3 data_format_trans.py new_10_genomes_tblastn.out new_10_genomes_tblastn.trans.txt s
```
```
python3 data_format_trans.py META23IMJSL03-10_tblastn_results_all.out META23IMJSL03-10_tblastn.trans.txt s
```

###Step 2, calculate the final results:
```
python3 gene_cluster_extract_tblastn.py gene_lable new_10_genomes_tblastn.trans.txt result/new_10_genomes_tblastn_output tblastn
```
```
python3 gene_cluster_extract_tblastn.py gene_lable META23IMJSL03-10_tblastn.trans.txt result/META23IMJSL03-10_tblastn_output tblastn
```