import pandas as pd
import numpy as np
import re
from Bio import Entrez
import time
import pickle
import os

print("🚀 开始处理数据...")

# =========================
# ⚠️ 必填
# =========================
Entrez.email = "nk3901@foxmail.com"

# --- 1. 读取数据 ---
df_rip = pd.read_csv('data/RIPs.csv')
df_ori = pd.read_csv('data/selected_ori_regions.csv', usecols=range(5))

print(f"RIP: {len(df_rip)}")
print(f"OriC: {len(df_ori)}")

# =========================
# 🧬 2. 先查 taxonomy（核心优化）
# =========================

from Bio import Entrez
import time

Entrez.email = "nk3901@foxmail.com"

def fetch_taxonomy(accessions, batch_size=200, max_retries=3):
    taxid_map = {}
    species_map = {}
    
    acc_list = list(accessions)
    
    for i in range(0, len(acc_list), batch_size):
        batch = acc_list[i:i+batch_size]
        
        for attempt in range(max_retries):
            try:
                handle = Entrez.esummary(
                    db="nuccore",
                    id=",".join(batch)
                )
                records = Entrez.read(handle)
                
                for rec in records:
                    acc = rec['AccessionVersion'].split('.')[0]
                    taxid = rec.get('TaxId', None)
                    title = rec.get('Title', '')
                    
                    taxid_map[acc] = taxid
                    species_map[acc] = title
                
                break  # 成功就跳出 retry
            
            except Exception as e:
                print(f"⚠️ batch {i} 第 {attempt+1} 次失败：{e}")
                time.sleep(2 + attempt * 2)  # 递增等待
        
        # 👉 每个 batch 之间休息（关键！）
        time.sleep(0.5)
        
        # 👉 进度显示
        if i % 2000 == 0:
            print(f"已完成 {i}/{len(acc_list)}")
    
    return taxid_map, species_map


print("🧬 获取 taxonomy（优化：只查 accession 一次）...")

# 只用 Ori 就够（推荐）
df_ori['Accession_Number'] = df_ori['Accession_Number'].str.split('.').str[0]
unique_acc = df_ori['Accession_Number'].dropna().unique()

print(f"需要查询 {len(unique_acc)} 个 accession")

# 缓存
cache_file = "taxonomy_cache.pkl"

if os.path.exists(cache_file):
    print("⚡ 使用缓存...")
    with open(cache_file, "rb") as f:
        taxid_map, species_map = pickle.load(f)
else:
    taxid_map, species_map = fetch_taxonomy(unique_acc)
    with open(cache_file, "wb") as f:
        pickle.dump((taxid_map, species_map), f)

# 映射到 Ori
df_ori['taxid'] = df_ori['Accession_Number'].map(taxid_map)
df_ori['species'] = df_ori['Accession_Number'].map(species_map)

# =========================
# 🌿 3. 物种清洗
# =========================

def clean_species(name):
    if pd.isna(name):
        return "unknown"
    
    name = name.strip()
    
    blacklist = [
        'uncultured', 'unidentified', 'metagenome',
        'environmental sample', 'bacterium', 'archaeon',
        'sp.', 'sp ', 'clone'
    ]
    
    if any(word in name.lower() for word in blacklist):
        return "unknown"
    
    name = re.sub(r' plasmid.*', '', name, flags=re.IGNORECASE)
    
    parts = name.split()
    
    if len(parts) == 1:
        return parts[0]
    
    genus, species = parts[0], parts[1]
    
    if not species.islower():
        return genus
    
    return f"{genus} {species}"

df_ori['species_clean'] = df_ori['species'].apply(clean_species)

# =========================
# 🧬 4. 构建 ID
# =========================

df_ori['ori_start'] = df_ori['Intergenic_Start']
df_ori['ori_end'] = df_ori['Intergenic_End']
df_ori['ori_id'] = (
    df_ori['Accession_Number'].astype(str) + "_" +
    df_ori['ori_start'].astype(str) + "_" +
    df_ori['ori_end'].astype(str)
)

df_rip['Accession_Number'] = df_rip['Accession_Number'].str.split('.').str[0]
df_rip['rep_start'] = df_rip['gene_start']
df_rip['rep_end'] = df_rip['gene_end']
df_rip['rep_id'] = (
    df_rip['Accession_Number'].astype(str) + "_" +
    df_rip['rep_start'].astype(str) + "_" +
    df_rip['rep_end'].astype(str)
)

# =========================
# 🔗 5. 配对（已经有 taxonomy）
# =========================

df_ori_small = df_ori[
    ['Accession_Number', 'ori_id', 'ori_start', 'ori_end', 'Intergenic_Sequence', 'species_clean', 'taxid']
]

df_rip_small = df_rip[
    ['Accession_Number', 'rep_id', 'rep_start', 'rep_end', 'gene_id', 'product', 'translation']
]

print("🔗 正在进行 Ori–Rep 配对...")

df_merged = pd.merge(
    df_ori_small,
    df_rip_small,
    on='Accession_Number',
    how='left'
)

print(f"配对后行数：{len(df_merged)}")

# =========================
# 📏 6. 距离
# =========================

ori_mid = (df_merged['ori_start'] + df_merged['ori_end']) / 2
rep_mid = (df_merged['rep_start'] + df_merged['rep_end']) / 2

df_merged['distance'] = abs(ori_mid - rep_mid)

# =========================
# 📦 7. 输出
# =========================

df_final = pd.DataFrame({
    'OriC sequence': df_merged['Intergenic_Sequence'],
    'ori_id': df_merged['ori_id'],
    'plasmid_id': df_merged['Accession_Number'],
    
    'species': df_merged['species_clean'],
    'taxid': df_merged['taxid'],
    
    'source': 'PLSDB',
    'pfamid_fast': np.nan,
    
    'rep_id': df_merged['rep_id'],
    'Rep_type_fast': df_merged['product'],
    'rep_seq': df_merged['translation'],
    
    'rep_dna_seq': np.nan,
    'full_replicon_seq': np.nan,
    'split': np.nan,
    '__index_level_0__': np.nan,
    
    'ori_start': df_merged['ori_start'],
    'ori_end': df_merged['ori_end'],
    'rep_start': df_merged['rep_start'],
    'rep_end': df_merged['rep_end'],
    'distance': df_merged['distance']
})

output_file = 'merged_final_optimized.csv'
df_final.to_csv(output_file, index=False)

print("✅ 完成！")
print(f"最终数据量：{len(df_final)}")
print(f"输出文件：{output_file}")

# =========================
# 📊 检查
# =========================
print("\n📊 检查：")
print("taxid 缺失：", df_final['taxid'].isna().mean())
print("unknown 物种比例：", (df_final['species'] == 'unknown').mean())