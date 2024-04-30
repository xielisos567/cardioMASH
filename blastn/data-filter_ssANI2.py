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

def dataFilter(data, setv, out):
	out_file = open(out, 'w')
	out_file.write('sample\tqseqid\tgene\tstrain\tspecies\ts_genome\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tgaps\tqcovs\tqcovhsp\tsseq\n')
	with open(data) as infile:
		for line in infile:
			line = re.sub('\[gene=.*?\].*\[gbkey=CDS\]_', '', line)
			line_info = line.strip().split('\t')
			identity = float(line_info[3])
			#print (line_info[1])
			try:
				gene_info, species_info = re.search('.*\[(.*?)\]_{0,}\[(.*)\]$', line_info[1]).groups()
			except:
				species_info = re.search('.*_\[(.*)\]_{0,}$', line_info[1]).group(1)
				try:
					gene_info = re.search('.*?_\d{1,}_\d{1,}_(.*?)_\[', line_info[1]).group(1)
				except:
					gene_info = "_".join(line_info[1].split('[')[0].split('_')[2:])[:-1]
				#species_info = re.search('.*_\[(.*)\]$', line_info[1]).group(1)
				#gene_info = "_".join(line_info[1].split('[')[0].split('_')[2:])[:-1]
				#if re.search('^[ctg|ODOSP|ALFI|GXM|HHO|CRH|PFJ|E0|NPD]', line_info[1]):
				#	gene_info = "_".join(line_info[1].split('[')[0].split('_')[2:])[:-1]
				#else:
				#	gene_info = ''
			if line_info[1].startswith('GCA'):
				genome_id = "_".join(line_info[1].split('_')[:5])
			elif re.search('^lcl', line_info[1]):
				genome_id = "_".join(line_info[1].split('_')[:6])
				#elif re.search('^[ERR|ctg|ALFI|ODOSP|GXM|HHO|CRH|NPD|BVU]', line_info[1]):
				#	genome_id = "_".join(line_info[1].split('_')[:2])
			else:
				genome_id = "_".join(line_info[1].split('_')[:2])
				#genome_id = ''
			species_info = re.sub("[\'\[\]]", "", species_info)
			species_id = "_".join(species_info.split('_')[:2])
			#species_id = re.sub("[\'\[\]]", "", species_id)
			#species_id = species_id.replace('\'','')
			if species_id not in setv:
				if identity >= 99:
					line_info.insert(2, species_id)
					line_info.insert(2, species_info)
					line_info.insert(2,gene_info)
					line_info.insert(2, genome_id)
					del line_info[1]
					out_data = '\t'.join(line_info)
					out_file.write(f'{out_data}\n')
			elif identity >= setv[species_id]:
				line_info.insert(2, species_id)
				line_info.insert(2, species_info)
				line_info.insert(2,gene_info)
				line_info.insert(2, genome_id)
				del line_info[1]
				out_data = '\t'.join(line_info)
				out_file.write(f'{out_data}\n')
	out_file.close()

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print ('Usage: python3 data-filter.py ssANI sample.out sample.screen.out')
		sys.exit()
	ssani, sampledata, sampleout = sys.argv[1:]
	setv = getSet(ssani)
	dataFilter(sampledata, setv, sampleout)
