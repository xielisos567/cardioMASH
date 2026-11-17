#step 0. makedb:
cd /pathway/to/genomes/
for i in `ls ./`
do
        makeblastdb -in ./$i -dbtype nucl
done
cd ..

#step 1. TBLASTN through the protein database of cvdMGC:
mkdir tblastn_result
for i in `less ./genomes.txt`
do
	tblastn -query cvdMGC/all_redundant_proteins.faa -db ./genomes/$i -evalue 1e-5 -qcov_hsp_perc 80 -max_target_seqs 1 -out tblastn_result/${i}_10.bls.fmt6 -outfmt "6 std gaps qcovs qcovhsp sseq" -num_threads 15
done

#step 2. Merge sequence alignment results:
cd tblastn_result
for i in `ls ./`; do cat ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../tblastn_genome-pep_results_all.out; done
cd ..

#step 3. Data format conversion:
mkdir final_results
python3 data_format_trans.py tblastn_genome-pep_results_all.out genomes_tblastn.trans.txt s

#step 4. Calculate the final results:
python3 gene_cluster_extract_tblastn.py gene_lable genomes_tblastn.trans.txt final_results/genomes_tblastn_output tblastn
