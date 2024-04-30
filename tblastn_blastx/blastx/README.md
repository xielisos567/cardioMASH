###blastx: sequence alignment of the nucleotide coding sequence to the eCAMGC protein database, the results are filtered by 40％ identity and 80％ coverage:
```
for i in `ls cds_seq/`
do
	time diamond blastx --db  db/all_CVD_translated_pep3.faa --query cds_seq/${i}_cds.fna --out result/${i}_CVD_40_80.tab --more-sensitive --evalue 1e-5 --id 40 --outfmt 6 --query-cover 80
done
```
###Step 1, data format conversion(letter c represents the format of the second column with both genome and gene numbers):
```
python3 data_format_trans2.py blastx_cds-pep_results_all.out blastx_cds-pep_results_all.trans.txt c
```
```
python3 data_format_trans2.py META23IMJSL03-10_BLASTX_all.out META23IMJSL03-10_BLASTX_all.trans.txt c
```

###Step 2, calculate the final results:
```
python3 gene_cluster_extract_blastx4.py gene_lable blastx_cds-pep_results_all.trans.txt result/blastx_cds-pep_results_all_output
```
```
python3 gene_cluster_extract_blastx4.py gene_lable META23IMJSL03-10_BLASTX_all.trans.txt result/META23IMJSL03-10_BLASTX_all_output
```