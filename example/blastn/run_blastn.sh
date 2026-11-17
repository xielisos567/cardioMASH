#step 0. makedb
cd /pathway/to/15_genomes
for i in `ls ./`
do
        makeblastdb -in ./$i -dbtype nucl
done
cd ..

#step1. BLASTN through the nucleotide database of cvdMGC
mkdir blastn_result
for i in `less ./15_genomes.txt`
do
        blastn -query ./cvdMGC/combine_CVD_cds.fasta -db ./15_genomes/${i}.fna -evalue 1e-5 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 1 -out blastn_result/${i}.bls.fmt6 -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 4
done

#step 2. Merge sequence alignment results
cd blastn_result
for i in `ls ./`; do cat ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../combine_blastn_results_all.out; done
cd ..

#step 3. ssANI filtering
mkdir final_result2
python3 data-filter_ssANI2.py ANI_result.txt combine_blastn_results_all.out 15_genomes_ssANI.txt

#step4. cardioMASH annotation
python3 gene_cluster_extract_final_ssANI6-fragment.py gene_lable 15_genomes_ssANI.txt final_result2/15_genomes_annotation
