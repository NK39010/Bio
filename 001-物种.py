import pandas as pd
from Bio import Entrez
import subprocess
import zipfile
import os
import re

# Step 1: 读取 CSV 并提取 plasmid_id 列
csv_file = "merged_final_optimized.csv"  # 填写你的 CSV 文件路径
df = pd.read_csv(csv_file)
plasmid_ids = df['plasmid_id'].dropna()  # 提取并去掉空值

# NCBI Entrez 设置
Entrez.email = "nk3901@foxmail.com"

# 用来存放结果的列表
results = []

# 正则表达式匹配 GCF 或 GCA 编号
def extract_assembly_accession(record):
    match = re.search(r"(GCF_\d+\.\d+|GCA_\d+\.\d+)", record)
    if match:
        return match.group(0)
    return None

# 解压并提取物种和复制方式
def extract_species_and_replication_type(zip_file):
    # 解压文件
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        zip_ref.extractall("temp_data")  # 解压到 temp_data 文件夹

    # 假设我们从解压后的文件中提取物种和复制方式
    # 这里可以根据解压后的文件内容来选择适当的文件进行解析
    species = "未知物种"  # 示例，需根据实际文件内容提取
    replication_type = "未知复制方式"  # 示例，需根据实际文件内容提取

    # 假设解压的文件是一个 XML 或 TXT 文件，获取其中的信息
    for file_name in os.listdir("temp_data"):
        if file_name.endswith(".xml"):  # 假设文件是 XML 格式
            with open(f"temp_data/{file_name}", "r") as file:
                content = file.read()
                # 解析文件内容以提取物种和复制方式
                # 示例：你可以根据 XML 文件的结构提取相关信息
                # 在这里，你需要编写具体的解析逻辑，可能是正则或 XML 解析器

                # 假设我们找到了相关信息
                species = "Example Species"
                replication_type = "Circular"  # 或者 "Linear"

    return species, replication_type

# Step 2: 获取组装号和其他信息
for plasmid_id in plasmid_ids:
    # Step 2.1: 通过 NC 编号查找组装号
    try:
        handle = Entrez.efetch(db="nucleotide", id=plasmid_id, rettype="gb", retmode="text")
        record = handle.read()
        
        # 提取组装号
        assembly_accession = extract_assembly_accession(record)
        if not assembly_accession:
            print(f"未找到组装号 for plasmid_id {plasmid_id}")
            continue  # 如果没有找到组装号，跳过当前 plasmid_id

        # Step 2.2: 用 datasets 调用获取物种和复制方式
        cmd = f"datasets download genome accession {assembly_accession} --filename temp_data.zip"
        subprocess.run(cmd, shell=True)

        # 解压并提取物种和复制方式
        species, replication_type = extract_species_and_replication_type("temp_data.zip")

        # Step 3: 存储结果
        results.append({"plasmid_id": plasmid_id, "assembly_accession": assembly_accession, "species": species, "replication_type": replication_type})
    
    except Exception as e:
        print(f"Error with plasmid_id {plasmid_id}: {e}")

# Step 4: 保存结果到 CSV
result_df = pd.DataFrame(results)
result_df.to_csv("output.csv", index=False)
print("结果已保存到 output.csv")