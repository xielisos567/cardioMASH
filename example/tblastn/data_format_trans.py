#Before using the python script, users should replace the space bar in the database with an underscore and then perform sequence alignment.
#!/bin/python3

import re
import os
import sys

def getSet(ssani):
	set_value = {}
	with open(ssani) as infile:
		next(infile)
		for line in infile:
			if re.search('^$', line):
				continue
			strain, value = line.strip().split('\t')
			set_value[strain.strip().replace(' ','_')] = float(value)
	return set_value

def dataFormatTrans(data, out, dtype):
	out_file = open(out, 'w')
	out_file.write('sample\tqseqid\tgene\tstrain\tspecies\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tgaps\tqcovs\tqcovhsp\tsseq\n')
	with open(data) as infile:
		for line in infile:
			line = re.sub('\[protein_id=.*?\].*\[gbkey=CDS\]_', '', line)
			line_info = line.strip().split('\t')
			identity = float(line_info[3])
			if dtype == 'c':
				gene_species = line_info[2].strip('_')
			else:
				gene_species = line_info[1].strip('_')
			#print (gene_species)
			try:
				gene_info, species_info = re.search('.*\[(.*?)\]_\[(.*)\]$', gene_species).groups()
			except:
				#print (line)
				species_info = re.search('.*_\[(.*)\]$', gene_species).group(1)
				gene_info = "_".join(gene_species.split('[')[0].split('_')[2:])[:-1]
				#if re.search('^[ctg|ODOSP|ALFI|GXM|HHO|CRH|PFJ|E0|NPD]', line_info[1]):
				#	gene_info = "_".join(line_info[1].split('[')[0].split('_')[2:])[:-1]
				#else:
				#	gene_info = ''
			if gene_species.startswith('GCA'):
				genome_id = "_".join(gene_species.split('_')[:5])
			elif re.search('^lcl', gene_species):
				genome_id = "_".join(gene_species.split('_')[:6])
				#elif re.search('^[ERR|ctg|ALFI|ODOSP|GXM|HHO|CRH|NPD|BVU]', line_info[1]):
				#	genome_id = "_".join(line_info[1].split('_')[:2])
			else:
				genome_id = "_".join(gene_species.split('_')[:2])
				#genome_id = ''
				if re.search('\d$', genome_id):
					pass
				elif re.search('\.\d{1,}' ,genome_id.split('_')[0]):
					pos = re.search('.*?(\d{1,})\.\d{1,}_', genome_id).group(1)
					genome_id = genome_id.split('_')[0]+"_"+pos
				else:
					gene_info = genome_id.split('_')[-1]+'_'+gene_info
					genome_id = '_'.join(genome_id.split('_')[:-1])
			if dtype == 'c':
				pass
			else:
				try:
					genome_id1 = re.search('(.*)\.fna\.bls\.', line_info[0]).group(1)
					genome_id = genome_id1+'_'+genome_id
				except:
					pass
			species_info = re.sub("[\"\'\[\]]", "", species_info)
			species_id = "_".join(species_info.split('_')[:2])
			#species_id = re.sub("[\'\[\]]", "", species_id)
			#species_id = species_id.replace('\'','')
			line_info.insert(2, species_id)
			line_info.insert(2, species_info)
			line_info.insert(2,gene_info)
			line_info.insert(2, genome_id)
			
			del line_info[1]
			#del gene_species
			out_data = '\t'.join(line_info)
			out_file.write(f'{out_data}\n')
	out_file.close()

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print ('Usage: python3 data-filter.py infile trans_file dtype')
		print ('\tdtype: s for caiTABCDE and c for others')
		sys.exit()
	indata, transdata, dtype = sys.argv[1:]
	dataFormatTrans(indata, transdata, dtype)
