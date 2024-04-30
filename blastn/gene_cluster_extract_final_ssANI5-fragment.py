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
	代码出现在rmDup函数中，对应代码位置在247-255行，修改对应gene_lable下的数值就可以
'''

# 获取基因id
def getGeneID(geneinfo, tables):
	if '/' in geneinfo:
		type_info = geneinfo.split('/')
		for ginfo in type_info:
			if ':' in ginfo:
				gtype, gene_id_tmp = ginfo.split(':')
				#gene_id_tmp = re.findall(':{0,}_{0,}([a-zA-Z_0-9\-]{1,})', ginfo)
				gene_id_tmp = gene_id_tmp.strip().strip('_')
				if gene_id_tmp in tables[gtype]:
					gene_id = gene_id_tmp
					gene_type = gtype
			else:
				gene_id = ginfo.strip()
				for gtype in tables:
					if gene_id in tables[gtype]:
						gene_type = gtype
						break
	else:
		if ':' in geneinfo:
			gtype, gene_id_tmp = geneinfo.split(':')
			#gene_id_tmp = re.findall(':{0,}_{0,}([a-zA-Z_0-9\-]{1,})', geneinfo)
			gene_id_tmp = gene_id_tmp.strip().strip('_')
			if gene_id_tmp in tables[gtype]:
				gene_id = gene_id_tmp
				gene_type = gtype
		else:
			gene_id = geneinfo.strip()
			for gtype in tables:
				if gene_id in tables[gtype]:
					gene_type = gtype
					break
	return gene_id, gene_type

# 定义函数，返回一个 defaultdict(list) 对象
def defaultdict_list():
	return defaultdict(list())

# 读取表格文件，将表格中的基因信息存储到字典 genes 中，并将表格信息存储到 defaultdict(list) 对象 tables 中
def getTable(file, name, table, genes, table_count):
	with open(file) as infile:
		for line in infile:
			extract_info = re.split(':|/', line)[-1].strip()
			table[name].append(extract_info)
			table_count[name] += 1
			genes[extract_info] = name
	return table, genes, table_count

# 读取数据文件，将符合条件的数据存储到输出文件 out 中
def getAllOk1(data, tables, genes, out):
	data_file = open(data)
	fail_file = open(f'{out}.fail', 'w')
	headers = next(data_file).strip().split('\t')
	out_fail_header = '\t'.join(headers)
	fail_file.write(out_fail_header+'\n')
	data_gene_info = defaultdict(list)
	data_line_num_info = defaultdict(list)
	line_num = 1
	# 将数据文件中的基因信息存储到字典 data_gene_info 中，将行号存储到字典 data_line_num_info 中
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
		#gene_id_tmp = re.split(':',re.split('/', line_info[2])[0])[-1].strip().strip('_')
		species_id = line_info[3].strip()
		key_id = ';'.join([sample_id, species_id, genome_id, type_lable])
		data_line_num_info[key_id].append(line_num)
		data_gene_info[key_id].append(gene_id)
		line_num += 1
	data_file.close()
	ok_line_list = []
	#print (data_gene_info['META23IMJSL03-18_bin_3.fasta.bls.fmt6;Escherichia_coli_2362-75_strain=2362-75;lcl|NZ_QTXO01000002.1_cds'])
	# 遍历 tables 中的表格信息，查找符合条件的数据行号
	for table in tables:
		for key in data_gene_info:
			data_gene_list = set(data_gene_info[key])
			table_list = set(tables[table])
			if len(table_list.intersection(data_gene_list)) == len(table_list):
				#print (table_list, data_gene_list)
				ok_line_list.extend(data_line_num_info[key])
			elif table in ['acetate2butyrate', 'glutamate2butyric'] and len(list(table_list.intersection(data_gene_list))) >= 5:
				ok_line_list.extend(data_line_num_info[key])
			elif table == 'p-cresol' and len(list(set(['hpdA', 'hpdB']).intersection(data_gene_list))) >= 2:
				ok_line_list.extend(data_line_num_info[key])
			elif table == 'LPS' and len(list(table_list.intersection(data_gene_list))) >= 6:
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
	# 将符合条件的数据写入输出文件 out 中
	for line in data_file:
		if line_num in ok_line_list:
			line_info = line.strip().split('\t')
			if line_info[2] == '':
				continue
			#gene_id = re.split(':|/', line_info[2])[-1].strip()
			#print (line, line_info[2])
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
	headers.insert(1, 'qseq_type')
	header = '\t'.join(headers)
	out_file = open(f'{out}.filter.xls', 'w')
	out_file.write(f'{header}\n')
	line_num = 1
	# 将符合条件的数据写入输出文件 out 中
	for line in data_file:
		line_info = line.strip().split('\t')
		if line_info[2] == '':
			fail_file.write(line)
			continue
		#gene_id = re.split(':|/', line_info[2])[-1].strip()
		#print (line, line_info[2])
		gene_id, type_lable = getGeneID(line_info[2], tables)
		line_info.insert(2, genes[gene_id])
		if line_info[1].startswith('lcl'):
			genome_id = '_'.join(line_info[1].split('_')[:-3]).strip()
		else:
			genome_id = '_'.join(line_info[1].split('_')[:-1]).strip()
		line_info.insert(1, genome_id)
		line = '\t'.join(line_info)
		out_file.write(f'{line}\n')
	out_file.close()
	data_file.close()
	fail_file.close()

# 从列表 data 中删除列表 lst 中的元素，并返回删除后的列表
def delLst(data, lst):
	for item in lst:
		# data.remove(item)
		# data = [x for x in data if x not in lst and not pd.isna(x)]
		data = [x for x in data if not x.equals(item) and not pd.isna(x).any().any()]
	return data

# 获取连续基因的数量和去重后的数据
def getConsCount(df, generange):
	df.drop_duplicates(inplace=True)
	cols = df.columns
	value_list = []
	del_value = []
	#total_count = 1
	# 获取基因编号列表 value_list，并将其排序
	for data in df['qseqid']:
		try:
			value = re.search('_[a-zA-Z]{0,}(\d{1,})$', data).group(1)
		except:
			value = re.search('(\d{1,})_.*$', data).group(1)
		value_list.append(int(value))
	value_list.sort()
	start = value_list[0]
	cons_count = 1
	# 计算连续基因的数量
	for item in value_list[1:]:
		#total_count += 1
		if item - start < generange:
			cons_count += 1
			start = item
		else:
			if value_list.index(start) == 0:
				del_value.append(start)
				start = item
			else:
				del_value.append(item)
			#start = item
	df_new = pd.DataFrame()
	# 将去重后的数据存储到 DataFrame 对象 df_new 中
	for index, row in df.iterrows():
		data = row['qseqid']
		try:
			value = re.search('_[a-zA-Z]{0,}(\d{1,})$', data).group(1)
		except:
			value = re.search('(\d{1,})_.*$', data).group(1)
		#value = re.search('_[a-zA-Z]{0,}(\d{1,})$', data).group(1)
		if int(value) in del_value:
			continue
		else:
			df_new = df_new.append(row)
	df_new = df_new[cols]
	#df_new.index = range(len(df_new))
	df_new.drop_duplicates(inplace=True)
	#print (len(df_new), len(df_new1))
	#return total_count, cons_count, df_new
	return cons_count, df_new

# 去除重复数据
def rmDup(out_body, table_count, generange=10):
	data = f'{out_body}.filter.xls'
	out = f'{out_body}.complete.rmdup.xls'
	in_data = pd.read_table(data, header=0)
	in_data = in_data.drop_duplicates(subset=['sample', 'qseq_type', 'Type', 'gene', 'sstart', 'send'])
	#print (len(in_data))
	#in_data.to_csv('test.xls', sep='\t', index=False)
	all_grouped_data = in_data.groupby(['sample', 'Type', 's_genome', 'qseq_type'])
	print (len(all_grouped_data))
	#all_grouped_data = in_data.groupby(['sample', 'Type'])
	rmdup_data = defaultdict(list)
	rmdup_data_tmp = defaultdict(list)
	fail_df = pd.DataFrame()
	# 遍历分组后的数据，去除重复数据
	total = 0
	for group_data in all_grouped_data:
		#print (group_data[1][['s_genome', 'gene']])
		total += 1
		if total % 10000 == 0:
			print (f'analysys {total}')
		gene_type = group_data[0][1]
		#data_id = ';'.join(group_data[0][:-1])
		data_id = ';'.join(group_data[0][:-1])
		#rmdup_group = group_data[1].drop_duplicates(subset=['sample', 'gene', 'Type', 'sstart', 'send', 's_genome'])
		group_count, group_df = getConsCount(group_data[1], generange)
		#if gene_type == 'p-cresol' and len(set(group_df['gene'].to_list()).intersection(set(['hpdA', 'hpdB']))) < 2: #修改p-cresol
		#	fail_df = fail_df.append(group_df)
		#	continue
		#if group_data[0][0] == 'HF07_350_bin_12.fasta.bls.fmt6' and group_data[0][2]=='ctg24' and group_count == 4:
		#	print (group_count, table_count[gene_type])
		#	print (group_df)
		#	group_df.to_csv('test.porA.xls', sep='\t')
		if gene_type not in ['acetate2butyrate', 'glutamate2butyric']:
			if gene_type == 'LPS':
				if group_count < 6: #修改LPS数量
					fail_df = fail_df.append(group_df)
					continue
			elif gene_type == 'p-cresol':
				if len(set(group_df['gene'].to_list()).intersection(set(['hpdA', 'hpdB']))) < 2:
					#print (gene_type, group_df, set(group_df['gene'].to_list()).intersection(set(['hpdA', 'hpdB'])))
					fail_df = fail_df.append(group_df)
					continue
			elif group_count < table_count[gene_type]:
				fail_df = fail_df.append(group_df)
				continue
		elif group_count < 5:
			fail_df = fail_df.append(group_df)
			continue
		if gene_type == 'p-cresol':
			print (group_df)
		if data_id not in rmdup_data:
			rmdup_data[data_id].append(group_df)
		else:
			#drop_list = []
			add = 0
			for candicate_data in rmdup_data_tmp[data_id]:
				#rm_data_del = defaultdict(list)
				start = sorted(group_df['sstart'])
				end = sorted(group_df['send'])
				instart = sorted(candicate_data['sstart'])
				inend = sorted(candicate_data['send'])
				# 判断两个数据是否重复
				if start == instart or end == inend:
					#total_count, cand_count, cand_df = getConsCount(candicate_data, generange)
					#cand_count, cand_df = getConsCount(candicate_data, generange)
					if candicate_data['mismatch'].sum() > group_df['mismatch'].sum():
						rmdup_data[data_id].remove(candicate_data)
						add = 1
						'''
						for item in rmdup_data[data_id]:
							if all(item[['sample', 'Type', 's_genome', 'sstart', 'send']]==cand_df[['sample', 'Type', 's_genome', 'sstart', 'send']]):
							#if all(item[['sstart', 'send']]==cand_df[['sstart', 'send']]):
								continue
							else:
								rm_data_del[data_id].append(item)
						#rmdup_data[data_id].remove(cand_df)
						rmdup_data[data_id] = rm_data_del[data_id]
						rmdup_data[data_id].append(group_df)
						'''
				elif set(start).intersection(set(instart)) == set(instart) or set(end).intersection(set(inend)) == set(inend):
					rmdup_data[data_id].remove(candicate_data)
					add = 1
			if add == 1:
				rmdup_data[data_id].append(group_df)
		rmdup_data_tmp = rmdup_data
		#rmdup_data[data_id] = delLst(rmdup_data[data_id], drop_list)
	out_data = pd.DataFrame()
	# 将去重后的数据存储到 DataFrame 对象 out_data 中，并将其写入输出文件 out 中
	for key in rmdup_data:
		for item in rmdup_data[key]:
			out_data = out_data.append(item)
	out_data.to_csv(out, sep='\t', index=False)
	fail_df.to_csv(f'{out_body}.fragmentary.xls', sep='\t', index=False)
	return out_data, fail_df

def rmDup1(data, out, table_count, generange=10):
	in_data = pd.read_table(data, header=0)
	all_grouped_data = in_data.groupby(['sample', 'Type', 'Genome'])
	print (len(all_grouped_data))
	#all_grouped_data = in_data.groupby(['sample', 'Type'])
	rmdup_data = defaultdict(list)
	rmdup_data_tmp = defaultdict(list)
	# 遍历分组后的数据，去除重复数据
	total = 0
	for group_data in all_grouped_data:
		total += 1
		if total % 100 == 0:
			#print (f'analysys {total}')
			pass
		gene_type = group_data[0][1]
		#data_id = ';'.join(group_data[0][:-1])
		rmdup_group = group_data[1].drop_duplicates(subset=['sample', 'Type', 'Genome', 'sstart', 'send'])
		data_id = ';'.join(group_data[0][:-2])
		group_count, group_df = getConsCount(rmdup_group, generange)
		if gene_type not in ['acetate2butyrate', 'glutamate2butyric']:
			if gene_type == 'LPS' and group_count < 6:
				continue
			elif gene_type == 'p-cresol' and group_count < 2:
				continue
			elif group_count < table_count[gene_type]:
				continue
		elif group_count < 5:
			continue
		if data_id not in rmdup_data:
			rmdup_data[data_id].append(group_df)
		else:
			#drop_list = []
			for candicate_data in rmdup_data_tmp[data_id]:
				rm_data_del = defaultdict(list)
				start = sorted(group_df['sstart'])
				end = sorted(group_df['send'])
				instart = sorted(candicate_data['sstart'])
				inend = sorted(candicate_data['send'])
				# 判断两个数据是否重复
				if start == instart or end == inend:
					#total_count, cand_count, cand_df = getConsCount(candicate_data, generange)
					cand_count, cand_df = getConsCount(candicate_data, generange)
					if group_count > cand_count:
						#drop_list[cand_df] = ''
						for item in rmdup_data[data_id]:
							#if all(item[['sample', 'Type', 'Genome', 'sstart', 'send']]==cand_df[['sample', 'Type', 'Genome', 'sstart', 'send']]):
							if all(item[['sample', 'Type', 'Genome', 'sstart', 'send']]==cand_df[['sample', 'Type', 'Genome', 'sstart', 'send']]):
								continue
							else:
								rm_data_del[data_id].append(item)
						rmdup_data[data_id] = rm_data_del[data_id]
						rmdup_data[data_id].append(group_df)
					elif group_count == cand_count:
						if cand_df['mismatch'].sum() > group_df['mismatch'].sum():
							for item in rmdup_data[data_id]:
								if all(item[['sample', 'Type', 'Genome', 'sstart', 'send']]==cand_df[['sample', 'Type', 'Genome', 'sstart', 'send']]):
									continue
								else:
									rm_data_del[data_id].append(item)
							#rmdup_data[data_id].remove(cand_df)
							rmdup_data[data_id] = rm_data_del[data_id]
							rmdup_data[data_id].append(group_df)
			rmdup_data_tmp = rmdup_data
			#rmdup_data[data_id] = delLst(rmdup_data[data_id], drop_list)
	out_data = pd.DataFrame()
	# 将去重后的数据存储到 DataFrame 对象 out_data 中，并将其写入输出文件 out 中
	for key in rmdup_data:
		for item in rmdup_data[key]:
			out_data = out_data.append(item)
	out_data.to_csv(out, sep='\t', index=False)
	return out_data

# 统计数据信息
def getStat(df, out, tablelist):
	out_stat = pd.DataFrame(columns=tablelist)
	all_grouped_data = df.groupby(['sample', 'Type', 's_genome'])
	# 统计每个样本中每种基因的数量
	for group_data in all_grouped_data:
		sample_id, gene_type = group_data[0][:2]
		'''
		try:
			if not np.isnan(out_stat.loc[sample_id, gene_type]) :
				out_stat.loc[sample_id, gene_type] += 1
		except:
			out_stat.loc[sample_id, gene_type] = 1
		'''
		try:
			out_stat.loc[sample_id, gene_type] += 1
		except:
			out_stat.loc[sample_id, gene_type] = 1
		if np.isnan(out_stat.loc[sample_id, gene_type]):
			out_stat.loc[sample_id, gene_type] = 1
	out_stat = out_stat.fillna(0)
	out_stat.to_excel(out)

def rmDupUncomplete(data, out):
	#data_uncomplete = pd.read_table(data, header=0)
	all_grouped_data = data.groupby(['sample', 'Type', 's_genome', 'strain'])
	rmdup_data = defaultdict(list)
	rmdup_data_tmp = defaultdict(list)
	# 遍历分组后的数据，去除重复数据
	total = 0
	for group_data in all_grouped_data:
		total += 1
		if total % 10000 == 0:
			print (f'analysys {total}')
		gene_type = group_data[0][1]
		#data_id = ';'.join(group_data[0][:-1])
		data_id = ';'.join(group_data[0][:-2])
		#print (data_id)
		#rmdup_group = group_data[1].drop_duplicates(subset=['sample', 'Type', 'Genome', 'sstart', 'send'])
		#group_count, group_df = getConsCount(group_data[1], generange)
		#group_df = group_data[1].drop_duplicates()
		group_df = group_data[1].drop_duplicates(subset=['sample', 'Type', 'sstart', 'send'])
		#print (group_df)
		if data_id not in rmdup_data:
			rmdup_data[data_id].append(group_df)
			rmdup_data_tmp = rmdup_data
		else:
			#drop_list = []
			for candicate_data in rmdup_data_tmp[data_id]:
				#cand_genome = candicate_data['Genome'].iloc[0]
				#group_genome = group_df['Genome'].iloc[0]
				rm_data_del = defaultdict(list)
				start = group_df['sstart'].sort_values()
				end = group_df['send'].sort_values()
				instart = candicate_data['sstart'].sort_values()
				inend = candicate_data['send'].sort_values()
				# 判断两个数据是否重复
				if len(start) == len(instart):
					#total_count, cand_count, cand_df = getConsCount(candicate_data, generange)
					if candicate_data['pident'].sum() < group_df['pident'].sum():
						for item in rmdup_data[data_id]:
							item_start = group_df['sstart'].sort_values()
							if all(start == item_start):
								continue
							else:
								rm_data_del[data_id].append(item)
						#rmdup_data[data_id].remove(cand_df)
						rmdup_data[data_id] = rm_data_del[data_id]
						rmdup_data[data_id].append(group_df)
				elif len(start)>len(instart) and len(instart[instart.isin(start)]) == len(instart):
					#len(group_df['sstart'][group_df['sstart'].apply(lambda x:x in candicate_data['sstart'])]) == len(candicate_data['sstart']):
					#len(start)>len(instart) and len(set(start).intersection(set(instart))) == len(instart):
					for item in rmdup_data[data_id]:
						item_start = sorted(group_df['sstart'])
						if len(start)>len(item_start) and len(set(start).intersection(set(item_start))) == len(item_start):
							continue
						else:
							rm_data_del[data_id].append(item)
					rmdup_data[data_id] = rm_data_del[data_id]
					rmdup_data[data_id].append(group_df)
				elif len(start)<len(instart) and len(start[start.isin(start)]) == len(start):
					break
				else:
					rmdup_data[data_id].append(group_df)
			rmdup_data_tmp = rmdup_data
			#rmdup_data[data_id] = delLst(rmdup_data[data_id], drop_list)
	out_data = pd.DataFrame()
	# 将去重后的数据存储到 DataFrame 对象 out_data 中，并将其写入输出文件 out 中
	for key in rmdup_data:
		for item in rmdup_data[key]:
			out_data = out_data.append(item)
	out_data.to_csv(out, sep='\t', index=False)
	return out_data

# 统计数据信息
def getUncompleteStat(df, out, tablelist):
	#df = pd.read_table(df, header=0)
	out_stat = pd.DataFrame(columns=tablelist)
	all_grouped_data = df.groupby(['sample', 'Type'])
	# 统计每个样本中每种基因的数量
	for group_data in all_grouped_data:
		sample_id, gene_type = group_data[0][:2]
		'''
		try:
			if not np.isnan(out_stat.loc[sample_id, gene_type]) :
				out_stat.loc[sample_id, gene_type] += 1
		except:
			out_stat.loc[sample_id, gene_type] = 1
		'''
		try:
			out_stat.loc[sample_id, gene_type] += 1
		except:
			out_stat.loc[sample_id, gene_type] = 1
		if np.isnan(out_stat.loc[sample_id, gene_type]):
			out_stat.loc[sample_id, gene_type] = 1
	out_stat = out_stat.fillna(0)
	out_stat.to_excel(out)

def geneStat(indata, tables, out):
	class1 = {}
	class2 = []
	for key in tables:
		for item in tables[key]:
			class1[item] = key
			class2.append(item)
	out_stat = pd.DataFrame(columns=class2)
	data = pd.read_table(indata, header=0)
	out_data = pd.DataFrame()
	all_grouped_data = data.groupby(['sample', 's_genome'])
	for group_data in all_grouped_data:
		#gene_type = group_data[0][1]
		#data_id = ';'.join(group_data[0][:-1])
		#data_id = ';'.join(group_data[0])
		#print (data_id)
		#rmdup_group = group_data[1].drop_duplicates(subset=['sample', 'Type', 'Genome', 'sstart', 'send'])
		#group_count, group_df = getConsCount(group_data[1], generange)
		#group_df = group_data[1].drop_duplicates()
		group_df = group_data[1].drop_duplicates(subset=['sample','sstart', 'send'])
		group_df = group_df.drop_duplicates(subset=['sample','sstart'])
		group_df = group_df.drop_duplicates(subset=['sample','send'])
		out_data = out_data.append(group_df)
	out_data.to_csv(f'{out}.gene.data.xls', sep='\t', index=False)
	rmdup_data = open(f'{out}.gene.data.xls')
	next(rmdup_data)
	for line in rmdup_data:
			line_info = line.strip().split('\t')
			#sample_id = line_info[0].split('.')[0]
			sample_id = re.search('(.*bin.\d{1,}).', line_info[0]).group(1)
			for genetype in class2:
				if re.search(genetype, line_info[4]):
					gene = genetype
			'''
			try:
					gene = line_info[4].split('/')[1]
			except:
					gene = line_info[4]
			'''
			try:
					out_stat.loc[sample_id, gene] += 1
			except:
					out_stat.loc[sample_id, gene] = 1
			if np.isnan(out_stat.loc[sample_id, gene]):
					out_stat.loc[sample_id, gene] = 1
	rmdup_data.close()
	out_stat = out_stat.fillna(0)
	out_stat.to_excel(f'{out}.gene.stat.xlsx')

def geneStat1(indata, tables, out):
	class1 = {}
	class2 = []
	for key in tables:
		for item in tables[key]:
			class1[item] = key
			class2.append(item)
	out_stat = pd.DataFrame(columns=class2)
	data = pd.read_table(indata, header=0)
	all_grouped_data = data.groupby(['sample','strain'])
	rmdup_data = defaultdict(list)
	rmdup_data_tmp = defaultdict(list)
	# 遍历分组后的数据，去除重复数据
	total = 0
	for group_data in all_grouped_data:
		total += 1
		if total % 10000 == 0:
			print (f'analysys {total}')
		gene_type = group_data[0][1]
		#data_id = ';'.join(group_data[0][:-1])
		data_id = ';'.join(group_data[0][:-1])
		#print (data_id)
		#rmdup_group = group_data[1].drop_duplicates(subset=['sample', 'Type', 'Genome', 'sstart', 'send'])
		#group_count, group_df = getConsCount(group_data[1], generange)
		#group_df = group_data[1].drop_duplicates()
		group_df = group_data[1].drop_duplicates(subset=['sample', 'sstart', 'send'])
		if data_id not in rmdup_data:
			rmdup_data[data_id].append(group_df)
			rmdup_data_tmp = rmdup_data
		else:
			#drop_list = []
			for candicate_data in rmdup_data_tmp[data_id]:
				#cand_genome = candicate_data['Genome'].iloc[0]
				#group_genome = group_df['Genome'].iloc[0]
				rm_data_del = defaultdict(list)
				start = sorted(group_df['sstart'])#.sort_values().reset_index()
				end = sorted(group_df['send'])#.sort_values().reset_index()
				instart = sorted(candicate_data['sstart'])#.sort_values().reset_index()
				inend = sorted(candicate_data['send'])#.sort_values().reset_index()
				# 判断两个数据是否重复
				if (start==instart or end==inend) and candicate_data['pident'].sum() < group_df['pident'].sum():
					#print (start, instart)
					for item in rmdup_data[data_id]:
						if len(item) == len(candicate_data):
							if all(item.reset_index() == candicate_data.reset_index()):
								continue
						else:
							rm_data_del[data_id].append(item)
					rmdup_data[data_id] = rm_data_del[data_id]
					rmdup_data[data_id].append(group_df)
					#total_count, cand_count, cand_df = getConsCount(candicate_data, generange)
					'''
					for item in rmdup_data[data_id]:
						item_start = sorted(group_df['sstart'])
						item_end = sorted(group_df['send'])
						if len(start)==len(item_start):
							if all(start==item_start) or all(end==item_end)):
							continue
						elif len(start)>len(item_start) and (len(set(start).intersection(set(item_start))) == len(item_start) or len(set(end).intersection(set(item_end))) == len(item_end)):
							continue
						elif len(start)<len(item_start) and (len(item_start[item_start.isin(start)]) == len(start) or len(item_end[item_end.isin(end)]) == len(end)):
							continue
						else:
							rm_data_del[data_id].append(item)
					rmdup_data[data_id] = rm_data_del[data_id]
					rmdup_data[data_id].append(group_df)
					'''
					#else:
					#	if candicate_data['pident'].sum() < group_df['pident'].sum():
					#		rm_data_del[data_id].append(item)
					#rmdup_data[data_id].remove(cand_df)
					#rmdup_data[data_id] = rm_data_del[data_id]
					#rmdup_data[data_id].append(group_df)
				elif len(start)>len(instart) and set(instart).intersection(set(start)) == set(instart) or set(inend).intersection(set(end)) == set(inend):
					#len(group_df['sstart'][group_df['sstart'].apply(lambda x:x in candicate_data['sstart'])]) == len(candicate_data['sstart']):
					#len(start)>len(instart) and len(set(start).intersection(set(instart))) == len(instart):
					for item in rmdup_data[data_id]:
						item_start = sorted(item['sstart'])
						item_end = sorted(item['send'])
						if len(start)>len(item_start) and (len(set(start).intersection(set(item_start))) == len(item_start) or len(set(end).intersection(set(item_end))) == len(item_end)):
							continue
						else:
							rm_data_del[data_id].append(item)
					rmdup_data[data_id] = rm_data_del[data_id]
					rmdup_data[data_id].append(group_df)
				elif len(start)<len(instart) and set(instart).intersection(set(start)) == set(start) or set(inend).intersection(set(end)) == set(end):
					break
				else:
					rmdup_data[data_id].append(group_df)
			rmdup_data_tmp = rmdup_data
			#rmdup_data[data_id] = delLst(rmdup_data[data_id], drop_list)
	out_data = pd.DataFrame()
	# 将去重后的数据存储到 DataFrame 对象 out_data 中，并将其写入输出文件 out 中
	for key in rmdup_data:
		for item in rmdup_data[key]:
			out_data = out_data.append(item)
	out_data.drop_duplicates(subset=['sample', 'Type', 'sstart', 'send'], inplace=True)
	out_data.to_csv(f'{out}.data', sep='\t', index=False)
	rmdup_data = open(f'{out}.data')
	next(rmdup_data)
	for line in rmdup_data:
			line_info = line.strip().split('\t')
			sample_id = line_info[0].split('.')[0]
			try:
					gene = line_info[2].split('/')[1]
			except:
					gene = line_info[2]
			try:
					out_stat.loc[sample_id, gene] += 1
			except:
					out_stat.loc[sample_id, gene] = 1
			if np.isnan(out_stat.loc[sample_id, gene]):
					out_stat.loc[sample_id, gene] = 1
	rmdup_data.close()
	out_stat = out_stat.fillna(0)
	out_stat.to_excel(out)
	#os.remove(f'{out}.data')

if __name__ == '__main__':
	table_path, data, out = sys.argv[1:]
	# 获取表格信息
	tables = defaultdict(list)
	table_count = defaultdict(int)
	all_type = []
	genes = {}
	for file in glob.glob(f'{table_path}/*txt'):
		name = os.path.basename(file).split('.')[0]
		all_type.append(name)
		tables, genes, table_count = getTable(file, name, tables, genes, table_count)
	# 过滤数据
	getAllOk(data, tables, genes, out)
	# 去除重复数据
	final_data, fail_data = rmDup(out, table_count, 10)
	# 统计基因簇信息
	getStat(final_data, f'{out}.complete.stat.xlsx', all_type)
	#geneStat(f'{out}.complete.rmdup.xls', tables, f'{out}.complete.gene.stat.xlsx')
	uncomplete_data = rmDupUncomplete(fail_data, f'{out}.fragmentary.rmdup.xls')
	getUncompleteStat(uncomplete_data, f'{out}.fragmentary.stat.xlsx', all_type)
	# 统计基因数量
	geneStat(f'{out}.filter.xls', tables, out)