import os
import re
import sys
import glob
import random
import numpy as np
import pandas as pd
from collections import defaultdict

# 获取基因id
def getGeneID(geneinfo, tables):
	#print (geneinfo)
	#global result11
	gene_type = None
	#if result11 == 2310:
		#print("2222222222222")
	if '/' in geneinfo:
		type_info = geneinfo.split('/')
		for ginfo in type_info:
			if ':' in ginfo:
				gtype, gene_id_tmp = ginfo.split(':')
				if re.search('degradaion_caiTABCDE', gtype):
					gtype = 'carnitine_degradaion_caiTABCDE'
				#gene_id_tmp = re.findall(':{0,}_{0,}([a-zA-Z_0-9\-]{1,})', ginfo)
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
			#gene_id_tmp = re.findall(':{0,}_{0,}([a-zA-Z_0-9\-]{1,})', geneinfo)
			#print (geneinfo, geneinfo.split(':'),tables[gtype])
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

	#result11 = 1+result11
	#if result11 == 2311:
		#print(gene_id, gene_type)
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

# 从列表 data 中删除列表 lst 中的元素，并返回删除后的列表
def delLst(data, lst):
	for item in lst:
		# data.remove(item)
		# data = [x for x in data if x not in lst and not pd.isna(x)]
		data = [x for x in data if not x.equals(item) and not pd.isna(x).any().any()]
	return data

# 获取连续基因的数量和去重后的数据
def getConsCount(df, generange):
	df.drop_duplicates(subset=['sample', 'Type', 'Group', 'strain', 'gene'], inplace=True)
	cols = df.columns
	value_list = []
	del_value = []
	cons_count = 1
	#total_count = 1
	# 获取基因编号列表 value_list，并将其排序
	for data in df['Genome']:
		#print (data)
		value = re.search('_[a-zA-Z]{0,}(\d{1,})$', data).group(1)
		value_list.append(int(value))
	value_list.sort()
	start = value_list[0]
	# 计算连续基因的数量
	for item in value_list[1:]:
		#total_count += 1
		if item - start <= generange:
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
		data = row['Genome']
		value = re.search('_[a-zA-Z]{0,}(\d{1,})$', data).group(1)
		if int(value) in del_value:
			continue
		else:
			#print("122222%s" % row)
			df_new = df_new.append(row)
	df_new = df_new[cols]
	#df_new.index = range(len(df_new))
	df_new.drop_duplicates(subset=['sample', 'Type', 'Group', 'strain', 'gene'], inplace=True)
	#print (len(df_new), len(df_new1))
	#return total_count, cons_count, df_new
	return cons_count, df_new

def rmdf(data, df1):
	data_rmdup = []
	for item in data:
		if all(df1==item):
			continue
		else:
			data_rmdup.append(item)
	return data_rmdup

def delCons(file, out):
	data = pd.read_table(file, sep='\t', header=0)
	all_grouped_data = data.groupby(['sample', 'Genome'])
	del_cons = pd.DataFrame()
	for group_data in all_grouped_data:
		df_cons = group_data[1].sort_values(by=['pident', 'mismatch'], ascending=[False, True]).reset_index()
		del_cons = del_cons.append(df_cons.iloc[0,:])
	del_cons.iloc[:,1:].to_csv(out, sep='\t', index=False)

def getAllOk(data, tables, genes, out):
	data_file = open(data)
	fail_file = open(f'{out}.fail.xls', 'w')
	headers = next(data_file).strip().split('\t')
	out_fail_header = '\t'.join(headers)
	fail_file.write(out_fail_header+'\n')
	headers.insert(2, 'Type')
	header = '\t'.join(headers)
	#headers.insert(1, 's_genome')
	#header = '\t'.join(headers)
	out_file = open(f'{out}.filter.xls', 'w')
	out_file.write(f'{header}\n')
	line_num = 1
	# 将符合条件的数据写入输出文件 out 中
	for line in data_file:
		line_info = line.strip().split('\t')
		if line_info[3] == '':
			fail_file.write(line)
			continue
		#gene_id = re.split(':|/', line_info[2])[-1].strip()
		#print (line, line_info[2])
		gene_id, type_lable = getGeneID(line_info[4], tables)
		if type_lable == None:
			continue
		line_info.insert(2, genes[gene_id])
		'''
		if line_info[1].startswith('lcl'):
			genome_id = '_'.join(line_info[1].split('_')[:-3]).strip()
		elif line_info[1].startswith('GCF'):
			genome_id = '_'.join(line_info[1].split('_')[:4]).strip()
		else:
			genome_id = '_'.join(line_info[1].split('_')[:-1]).strip()
		'''
		#line_info.insert(1, line_info[2])
		line = '\t'.join(line_info)
		out_file.write(f'{line}\n')
	out_file.close()
	data_file.close()
	fail_file.close()

def groupData(file, out, generange=10):
	data = pd.read_table(file, header=0)
	#print (data['Genome'])
	data['s1'] = data['Genome'].apply(lambda x:x.split('_')[-1]).astype(int)
	data = data.sort_values(by=['sample', 'Group','s1'])
	data['s2'] = data['s1']-data['s1'].shift(1)
	#print (data[['Genome', 's1', 's2']])
	data.to_csv(f'{out}.shift.xls', sep='\t', index=False)
	indata = open(f'{out}.shift.xls')
	outdata = open(f'{out}.shift.group.xls', 'w')
	header = next(indata).strip().split('\t')
	type_index = header.index('Type')
	sample_index = header.index('sample')
	group_index = header.index('Group')
	s2_index = header.index('s2')
	header.insert(4, 'Rank')
	#out_header = '\t'.join(header[:-2])
	out_header = '\t'.join(header)
	outdata.write(out_header+'\n')
	type_tmp = ''
	sample_tmp = ''
	group_tmp = ''
	rank = 0
	for line in indata:
		line_info = line.strip().split('\t')
		type_info = line_info[type_index]
		sample_info = line_info[sample_index]
		group_info = line_info[group_index]
		if type_tmp==type_info and sample_info==sample_tmp and group_tmp==group_info:
			if float(line_info[s2_index]) <= generange:
				line_info.insert(4, 's'+str(rank))
			else:
				rank += 1
				line_info.insert(4, 's'+str(rank))
				type_tmp = type_info
				sample_tmp = sample_info
				group_tmp = group_info
		else:
			rank += 1
			line_info.insert(4, 's'+str(rank))
			type_tmp = type_info
			sample_tmp = sample_info
			group_tmp = group_info
		#outdata.write('\t'.join(line_info[:-2])+'\n')
		outdata.write('\t'.join(line_info)+'\n')
	outdata.close()
	return pd.read_table(f'{out}.shift.group.xls', header=0)

# 去除重复数据
def rmDup(data, out, table_count, generange=10):
	#in_data = pd.read_table(data, header=0)
	in_data = groupData(data, 'add_rank')
	all_grouped_data = in_data.groupby(['sample', 'Type', 'Rank'])
	rmdup_data = defaultdict(list)
	rmdup_data_tmp = defaultdict(list)
	# 遍历分组后的数据，去除重复数据
	#total = 0
	fail_df = pd.DataFrame()
	for group_data in all_grouped_data:
		gene_type = group_data[0][1]
		data_id = ';'.join(group_data[0])
		group_count, group_df = getConsCount(group_data[1], generange)
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
		if data_id not in rmdup_data:
			rmdup_data[data_id].append(group_df)
			rmdup_data_tmp = rmdup_data
		else:
			add = 0
			for candicate_data in rmdup_data_tmp[data_id]:
				if len(group_df)>len(candicate_data):
					rmdup_data[data_id] = rmdf(rmdup_data[data_id], candicate_data)
					add = 1
				elif len(group_df)<len(candicate_data):
					continue
				else:
					if all(candicate_data['Genome'].reset_index()==group_df['Genome'].reset_index()):
						if candicate_data['mismatch'].sum() > group_df['mismatch'].sum():
							rmdup_data[data_id] = rmdf(rmdup_data[data_id], candicate_data)
							add = 1
					else:
						add = 1
			if add == 1:
				rmdup_data[data_id].append(group_df)
		rmdup_data_tmp = rmdup_data
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
	all_grouped_data = df.groupby(['sample', 'Type', 'strain'])
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

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print (f'python3 {sys.argv[0]} table_path input out_body')
		sys.exit()
	table_path, input_file, out_body = sys.argv[1:]
	delCons(input_file, f'{out_body}.delcons.xls')
	tables = defaultdict(list)
	table_count = defaultdict(int)
	all_type = []
	genes = {}
	for file in glob.glob(f'{table_path}/*txt'):
		name = os.path.basename(file).split('.')[0]
		all_type.append(name)
		tables, genes, table_count = getTable(file, name, tables, genes, table_count)
	# 过滤数据
	getAllOk(f'{out_body}.delcons.xls', tables, genes, f'{out_body}')
	final_data = rmDup(f'{out_body}.filter.xls', f'{out_body}.rmdup.xls', table_count, generange=10)
	# 统计基因簇信息
	getStat(final_data, f'{out_body}.complete.stat.xlsx', all_type)