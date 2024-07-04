# -*- coding: UTF-8 -*-
# FileName  : refine_domain_seq.py
# Time      : 2024/7/3
# Author    : xian

import sys
from Bio import SeqIO

# 打印参数使用说明
print("\n\tparameter assignment\n\tpython script.py idfile infastafile outfastafile\n")

# 从命令行参数中获取输入和输出文件名
idfile = sys.argv[1]
infastafile = sys.argv[2]
outfastafile = sys.argv[3]

# 读取ID文件并存储ID到列表中
id_list = []
with open(idfile) as f:
    for line in f:
        id_list.append(line.strip())

# 解析输入的FASTA文件
in_sequences = SeqIO.parse(infastafile, "fasta")
# 打开输出的FASTA文件准备写入
out_sequences = open(outfastafile, "w")

# 处理每个序列对象
for seq_record in in_sequences:
    seq_id = seq_record.id
    for id in id_list:
        if id in seq_id:
            # 更新序列的显示ID
            seq_record.id = id
            seq_record.description = ""
            # 将符合条件的序列写入输出文件
            SeqIO.write(seq_record, out_sequences, "fasta")
            break

# 关闭输出文件
out_sequences.close()

print("\tID sequence retrieved\n\n")