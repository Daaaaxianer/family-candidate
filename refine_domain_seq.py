# -*- coding: UTF-8 -*-
# FileName  : refine_domain_seq.py
# Time      : 2024/7/3
# Author    : xian

import sys
from Bio import SeqIO

# 检查命令行参数的数量是否为4，否则打印用法并退出
if len(sys.argv) != 5:
    print("Usage: python script.py <hmmout> <protein.fasta> <outfile> <e-value>")
    sys.exit(1)

# 从命令行参数中获取输入文件和参数
hmmout_file = sys.argv[1]
protein_fasta = sys.argv[2]
outfile = sys.argv[3]
e_value_threshold = float(sys.argv[4])

# 解析蛋白质序列文件
hmmout = SeqIO.parse(protein_fasta, "fasta")
# 准备输出文件，以FASTA格式写入
domfasta = open(outfile, "w")

# 初始化一个字典，用于保存满足条件的序列
keep = {}

# 打开hmmout文件并读取
with open(hmmout_file) as f:
    for line in f:
        line = line.strip()
        # 跳过以#开头的注释行
        if line.startswith("#"):
            continue
        
        # 按空白字符分割行内容
        arr = line.split()
        # 如果e-value大于阈值，跳过该行
        if float(arr[6]) > e_value_threshold:
            continue
        
        target = arr[0]
        brr = (int(arr[17]), int(arr[18]))
        
        # 如果target不在字典中，则添加
        if target not in keep:
            keep[target] = brr

# 重新打开蛋白质序列文件，处理每个序列
with open(protein_fasta) as f:
    for seq_record in SeqIO.parse(f, "fasta"):
        id = seq_record.id
        # 如果序列ID在字典中，处理该序列
        if id in keep:
            start, end = keep[id]
            print(f"{id} {start} {end}")
            
            # 提取子序列，BioPython使用0-based索引
            new_seq = seq_record[start-1:end]
            new_seq.id = seq_record.id
            new_seq.description = seq_record.description
            
            # 将新的子序列写入输出文件
            SeqIO.write(new_seq, domfasta, "fasta")

# 关闭输出文件
domfasta.close()