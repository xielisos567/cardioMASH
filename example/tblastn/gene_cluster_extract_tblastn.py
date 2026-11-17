import os
import re
import sys
import glob
import random
import numpy as np
import pandas as pd
from collections import defaultdict

'''
基因个数修改：
	代码出现在rmDup函数中，对应代码位置在273-282行，修改对应gene_lable下的数值就可以
'''

# 获取基因id
def getGeneID(geneinfo, tables):
	geneinfo = re.sub('^\d{1,}_', '', geneinfo)
	gene_type = None
	if '/' in geneinfo:
		type_info = geneinfo.split('/')
		for ginfo in type_info:
			if ':' in ginfo:
				gtype, gene_id_tmp = ginfo.split(':')
				if re.search('degradaion_caiTABCDE', gtype):
					gtype = 'carnitine_degradaion_caiTABCDE'
				gene_id_tmp = gene_id_tmp.strip().strip('_')
				if gene_id_tmp in tables[gtype]:
					gene_id = gene_id_tmp
					gene_type = gtype
			else:
				gene_id = ginfo.strip()
				for gtype in tables:
					if re.search('degradaion_caiTABCDE', gtype):
						gtype = 'carnitine_degradaion_caiTABCDE'
					if gene_id in tables[gtype]:
						gene_type = gtype
						break
	else:
		if ':' in geneinfo:
			gtype, gene_id_tmp = geneinfo.split(':')
			if re.search('degradaion_caiTABCDE', gtype):
				gtype = 'carnitine_degradaion_caiTABCDE'
			gene_id_tmp = gene_id_tmp.strip().strip('_')
			if gene_id_tmp in tables[gtype]:
				gene_id = gene_id_tmp
				gene_type = gtype
		else:
			gene_id = geneinfo.strip()
			for gtype in tables:
				if gene_id in tables[gtype]:
					if re.search('degradaion_caiTABCDE', gtype):
						gtype = 'carnitine_degradaion_caiTABCDE'
					gene_type = gtype
					break
	return gene_id, gene_type

def defaultdict_list():
	return defaultdict(list())

def getTable(file, name, table, genes, table_count):
	with open(file) as infile:
		for line in infile:
			extract_info = re.split(':|/', line)[-1].strip()
			table[name].add(extract_info)
			table_count[name] += 1
			genes[extract_info] = name
	return table, genes, table_count

def getAllOk1(data, tables, genes, out):
	data_file = open(data)
	fail_file = open(f'{out}.fail', 'w')
	headers = next(data_file).strip().split('\t')
	out_fail_header = '\t'.join(headers)
	fail_file.write(out_fail_header+'\n')
	data_gene_info = defaultdict(list)
	data_line_num_info = defaultdict(list)
	line_num = 1
	for line in data_file:
		line_info = line.strip().split('\t')
		if line_info[2] == '':
			fail_file.write(line)
			continue
		sample_id = line_info[0].strip()
		if line_info[1].startswith('lcl'):
			genome_id = '_'.join(line_info[1].split('_')[:-3]).strip()
		else:
			genome_id = '_'.join(line_info[1].split('_')[:-1]).strip()
		gene_id, type_lable = getGeneID(line_info[2], tables)
		species_id = line_info[3].strip()
		key_id = ';'.join([sample_id, species_id, genome_id, type_lable])
		data_line_num_info[key_id].append(line_num)
		data_gene_info[key_id].append(gene_id)
		line_num += 1
	data_file.close()
	ok_line_list = []
	table_sets = {t: set(v) for t, v in tables.items()}
	data_sets = {k: set(v) for k, v in data_gene_info.items()}
	for table, tset in table_sets.items():
		for key, dset in data_sets.items():
			inter = tset & dset
			if len(inter) == len(tset):
				ok_line_list.extend(data_line_num_info[key])
			elif table in ['acetate2butyrate', 'glutamate2butyric'] and len(inter) >= 5:
				ok_line_list.extend(data_line_num_info[key])
			elif table == 'p-cresol' and len(inter) >= 2:
				ok_line_list.extend(data_line_num_info[key])
			elif table == 'LPS' and len(inter) >= 6:
				ok_line_list.extend(data_line_num_info[key])
	ok_line_list = set(ok_line_list)
	data_file = open(data)
	next(data_file)
	out_file = open(out, 'w')
	headers.insert(2, 'Type')
	header = '\t'.join(headers)
	headers.insert(1, 'Genome')
	header = '\t'.join(headers)
	out_file.write(f'{header}\n')
	line_num = 1
	for line in data_file:
		if line_num in ok_line_list:
			line_info = line.strip().split('\t')
			if line_info[2] == '':
				continue
			gene_id, type_lable = getGeneID(line_info[2], tables)
			line_info.insert(2, genes[gene_id])
			if line_info[1].startswith('lcl'):
				genome_id = '_'.join(line_info[1].split('_')[:-3]).strip()
			else:
				genome_id = '_'.join(line_info[1].split('_')[:-1]).strip()
			line_info.insert(1, genome_id)
			line = '\t'.join(line_info)
			out_file.write(f'{line}\n')
		line_num += 1
	out_file.close()
	data_file.close()

def getAllOk(data, tables, genes, out):
	data_file = open(data)
	fail_file = open(f'{out}.fail.xls', 'w')
	headers = next(data_file).strip().split('\t')
	out_fail_header = '\t'.join(headers)
	fail_file.write(out_fail_header+'\n')
	headers.insert(2, 'Type')
	header = '\t'.join(headers)
	headers.insert(1, 'Genome')
	header = '\t'.join(headers)
	out_file = open(f'{out}.filter.xls', 'w')
	out_file.write(f'{header}\n')
	line_num = 1
	for line in data_file:
		line_info = line.strip().split('\t')
		if line_info[2] == '':
			fail_file.write(line)
			continue
		gene_id, type_lable = getGeneID(line_info[2], tables)
		if type_lable == None:
			continue
		line_info.insert(2, genes[gene_id])
		if line_info[1].startswith('lcl'):
			genome_id = '_'.join(line_info[1].split('_')[:-3]).strip()
		elif line_info[1].startswith('GCF'):
			genome_id = '_'.join(line_info[1].split('_')[:4]).strip()
		else:
			genome_id = '_'.join(line_info[1].split('_')[:-1]).strip()
		line_info.insert(1, genome_id)
		line = '\t'.join(line_info)
		out_file.write(f'{line}\n')
	out_file.close()
	data_file.close()
	fail_file.close()

def delLst(data, lst):
	for item in lst:
		data = [x for x in data if not x.equals(item) and not pd.isna(x).any().any()]
	return data

def getConsCount(df, generange):
	df.drop_duplicates(subset=['sample', 'Type', 'sseqid', 'strain', 'qseqid'], inplace=True)
	df.drop_duplicates(subset=['sample', 'Type', 'sseqid', 'strain', 'gene'], inplace=True)
	cols = list(df.columns)
	nums = df['qseqid'].str.extract(r'_([a-zA-Z]*)(\d+)$')[1].astype(int)
	value_list = sorted(nums.tolist())
	del_value = []
	cons_count = 1
	if value_list:
		start = value_list[0]
		for item in value_list[1:]:
			if item - start <= generange:
				cons_count += 1
				start = item
			else:
				if start == value_list[0]:
					del_value.append(start)
					start = item
				else:
					del_value.append(item)
	mask = ~nums.isin(del_value)
	df_new = df[mask][cols].copy()
	df_new.drop_duplicates(subset=['sample', 'Type', 'sseqid', 'strain', 'qseqid'], inplace=True)
	df_new.drop_duplicates(subset=['sample', 'Type', 'sseqid', 'strain', 'gene'], inplace=True)
	return cons_count, df_new

def rmDup(data, out, table_count, blasttype, generange=10):
	in_data = pd.read_table(data, header=0)
	if blasttype == 'tblastn':
		all_grouped_data = in_data.groupby(['sample', 'Type', 'sseqid', 'Genome', 'strain'])
	else:
		all_grouped_data = in_data.groupby(['sample', 'Type', 'sseqid', 'Genome', 'strain'])
	rmdup_data = defaultdict(list)
	rmdup_data_tmp = defaultdict(list)
	total = 0
	group_index = {}
	for group_data in all_grouped_data:
		total += 1
		gene_type = group_data[0][1]
		data_id = ';'.join(group_data[0][:-2])
		group_count, group_df = getConsCount(group_data[1], generange)
		group_index[f'group{total}'] = group_df
		if gene_type == 'p-cresol' and len(set(group_df['gene'].to_list()).intersection(set(['hpdA', 'hpdB']))) < 2:
			continue
		if gene_type not in ['acetate2butyrate', 'glutamate2butyric']:
			if gene_type == 'LPS' and group_count < 6:
				continue
			elif group_count<table_count[gene_type] and gene_type!='LPS' and gene_type!='p-cresol':
				continue
		elif group_count < 5:
			continue
		if data_id not in rmdup_data:
			rmdup_data[data_id].append(f'group{total}')
		else:
			add = 0
			drop_list = []
			count = 0
			for candicate_id in rmdup_data_tmp[data_id]:
				candicate_data = group_index[candicate_id]
				count += 1
				rm_data_del = defaultdict(list)
				start = sorted(group_df['sstart'])
				end = sorted(group_df['send'])
				instart = sorted(candicate_data['sstart'])
				inend = sorted(candicate_data['send'])
				if start == instart or end == inend:
					if (candicate_data['mismatch'].sum()>group_df['mismatch'].sum()):
						rmdup_data[data_id].remove(candicate_id)
						add = 1
					else:
						add = 0
						break
				elif (abs(candicate_data['sstart'].sort_values()-group_df['sstart'].sort_values()).sum()<=50) or (abs(candicate_data['send'].sort_values()-group_df['send'].sort_values()).sum()<=50):
					add = 0
					break
				elif set(start).intersection(set(instart)) == set(instart) or set(end).intersection(set(inend)) == set(inend):
					rmdup_data[data_id].remove(candicate_id)
					add = 1
				elif set(instart).intersection(set(start)) == set(start) or set(inend).intersection(set(end)) == set(end):
					add = 0
					break
				else:
					drop_list.append('1')
			if add == 1 or (add==0 and len(drop_list)==count):
				rmdup_data[data_id].append(f'group{total}')
		rmdup_data_tmp = rmdup_data
	frames = []
	for key in rmdup_data:
		for item in rmdup_data[key]:
			frames.append(group_index[item])
	if frames:
		out_data = pd.concat(frames, ignore_index=True)
	else:
		out_data = pd.DataFrame()
	out_data.to_csv(out, sep='\t', index=False)
	return out_data

def getStat(df, out, tablelist, blasttype):
	if blasttype == 'tblastn':
		out_stat = (df.drop_duplicates(subset=['sample', 'Type', 'sseqid', 'strain'])
		              .groupby(['sample', 'Type'])
		              .size()
		              .unstack(fill_value=0))
	else:
		out_stat = (df.drop_duplicates(subset=['sample', 'Type', 'sseqid', 'strain'])
		              .groupby(['sample', 'Type'])
		              .size()
		              .unstack(fill_value=0))
	out_stat = out_stat.reindex(columns=tablelist, fill_value=0)
	out_stat.to_excel(out)

if __name__ == '__main__':
	if len(sys.argv) != 5:
		print (f'Usage:\n   python3 {sys.argv[0]} table_path infile outprefix blast_type')
		sys.exit()
	table_path, data, out, blasttype = sys.argv[1:]
	out_path = '/'.join(os.path.abspath(out).split('/')[:-1])
	if not os.path.exists(out_path):
		os.mkdir(out_path)
	tables = defaultdict(set)
	table_count = defaultdict(int)
	all_type = []
	genes = {}
	for file in glob.glob(f'{table_path}/*txt'):
		name = os.path.basename(file).split('.')[0]
		all_type.append(name)
		tables, genes, table_count = getTable(file, name, tables, genes, table_count)
	getAllOk(data, tables, genes, out)
	final_data = rmDup(f'{out}.filter.xls', f'{out}.rmdup.xls', table_count, blasttype, 10)
	getStat(final_data, f'{out}.stat.xlsx', all_type, blasttype)
