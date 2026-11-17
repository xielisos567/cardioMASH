#step 0. makedb：
cd /pathway/to/cvdMGC
diamond makedb --in all_redundant_proteins.faa --db all_redundant_proteins --threads 8
cd ..

#step 1. BLASTX through the protein database：
mkdir blastx_result
for i in `ls cds/`
do
	time diamond blastx --db cvdMGC/all_redundant_proteins --query cds/$i --out blastx_result/${i}_CVD_40_80.tab --more-sensitive --evalue 1e-5 --id 40 --outfmt 6 --query-cover 80
done

#step 2. Merge sequence alignment results：
cd blastx_result
for i in `ls ./`; do cat ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../blastx_cds-pep_results_all.out; done
cd ..

#step 3. Data format conversion:
mkdir final_results
python3 data_format_trans.py blastx_cds-pep_results_all.out blastx_cds-pep_results_all.trans.txt c

#step 4. Calculate the final results
python3 gene_cluster_extract_blastx.py gene_lable blastx_cds-pep_results_all.trans.txt final_results/blastx_cds-pep_results_all_output
