###Step 1, mkblastdb of the genomes:
```
for i in `ls ./`
do
	makeblastdb -in ./$i -dbtype nucl
done
```

###Step 2, blastn through the nucleotide database of cvdMGC:
```
for i in `ls ./`
do
	blastn -query all_cds.fna -db ./$i -evalue 1e-5 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 1 -out result/$i.bls.fmt6 -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 4
done
```

###Step 3, filter the alignment results through the ssANI threshold (decompress 20_genomes_results_all.out.tar.gz firstly):
```
python3 data-filter_ssANI2.py ANI_result.txt 20_genomes_results_all.out 20_genomes_ssANI_output.txt
```

###Step 4, architecture search and sequence extraction of the physically clustered genes(the CVD associated MGCs are listed in the gene_lable folder):
```
python3 gene_cluster_extract_final_ssANI5-fragment.py gene_lable 20_genomes_ssANI_output.txt result/final_20_genomes_ssANI_output_statics
```

###Attentions:
pandas version < 1.4.0
