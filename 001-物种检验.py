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
Entrez.tool = "script"

# --- 1. 读取数据 ---
df_rip = pd.read_csv('all_RIP.csv')
df_ori = pd.read_csv('all_selected_ori_regions.csv', usecols=range(5))

print(f"RIP: {len(df_rip)}")
print(f"OriC: {len(df_ori)}")

# =========================
# 🧬 2. 基础查询 taxonomy（原有逻辑）
# =========================
def fetch_taxonomy(accessions, batch_size=200, max_retries=5):
    taxid_map = {}
    species_map = {}
    failed_acc = []
    acc_list = list(accessions)
    
    for i in range(0, len(acc_list), batch_size):
        batch = acc_list[i:i+batch_size]
        for attempt in range(max_retries):
            try:
                handle = Entrez.esummary(db="nuccore", id=",".join(batch))
                records = Entrez.read(handle)
                for rec in records:
                    acc = rec['AccessionVersion'].split('.')[0]
                    taxid = rec.get('TaxId', None)
                    title = rec.get('Title', '')
                    taxid_map[acc] = taxid
                    species_map[acc] = title
                break
            except Exception as e:
                print(f"⚠️ batch {i} 第 {attempt+1} 次失败: {e}")
                time.sleep(3 + attempt * 2)
        else:
            failed_acc.extend(batch)
        time.sleep(0.6)
        if i % 2000 == 0:
            print(f"已完成 {i}/{len(acc_list)}")
    return taxid_map, species_map, failed_acc

# 处理 accession
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
    taxid_map, species_map, failed_acc = fetch_taxonomy(unique_acc)
    with open(cache_file, "wb") as f:
        pickle.dump((taxid_map, species_map), f)

df_ori['taxid'] = df_ori['Accession_Number'].map(taxid_map)
df_ori['species'] = df_ori['Accession_Number'].map(species_map)

# =========================
# 🌟 3. 【彻底修复】物种清洗函数（适配NCBI ORGANISM）
# =========================
def clean_species(name):
    """
    专门适配 NCBI ORGANISM 字段的物种名清洗函数
    输入：原始物种名（如 "Marinovum algicola", "Escherichia coli str. K-12", "Bacillus sp.")
    输出：(标准化双名法名称, 是否异常)
    规范：Genus species（属名首字母大写，种名全小写）
    """
    if pd.isna(name) or str(name).strip() == "":
        return "unknown", True
    
    name = str(name).strip()
    
    # 1. 黑名单：只过滤真正无效的条目，不影响正常物种
    blacklist = [
        'uncultured', 'unidentified', 'metagenome',
        'environmental sample', 'bacterium', 'archaeon',
        'clone', 'vector', 'unclassified', 'synthetic',
        'construct', 'phage', 'virus', 'plasmid', 'chromosome'
    ]
    if any(word in name.lower() for word in blacklist):
        return "unknown", True
    
    # 2. 彻底清理所有冗余后缀（strain, plasmid, genome等）
    name = re.sub(
        r'\s*(strain|subsp\.|subspecies|serovar|serotype|biovar|pathovar|genomovar|pv\.|sv\.|bv\.|str\.).*$',
        '',
        name,
        flags=re.IGNORECASE
    )
    name = re.sub(
        r'\s*(plasmid|complete genome|whole genome|chromosome|isolate|clone).*$',
        '',
        name,
        flags=re.IGNORECASE
    )
    name = re.sub(r'\s+', ' ', name).strip()  # 合并多余空格
    
    # 3. 分割核心名称
    parts = name.split()
    if len(parts) < 1:
        return "unknown", True
    
    # 4. 提取属名和种名，严格遵循双名法规范
    genus = parts[0].capitalize()  # 属名首字母大写
    is_abnormal = False
    
    # 情况1：只有属名（如 "Marinovum"）→ 标记异常
    if len(parts) == 1:
        return genus, True
    
    # 情况2：双名法（如 "Marinovum algicola"）→ 标准化
    species_part = parts[1].lower()  # 种名全小写
    
    # 5. 种名有效性校验
    invalid_species = ['sp', 'sp.', 'spp', 'spp.', 'cf.', 'aff.', 'nov.', 'spnov']
    if species_part in invalid_species or len(species_part) < 2 or not species_part.isalpha():
        is_abnormal = True
        return genus, is_abnormal
    
    # 6. 最终标准化：Genus species
    clean_name = f"{genus} {species_part}"
    return clean_name, is_abnormal

df_ori[['species_clean', 'is_species_abnormal']] = df_ori['species'].apply(
    lambda x: pd.Series(clean_species(x))
)

# =========================
# 🔍 4. 【核心新增】unknown → 查 Assembly → 取 ORGANISM
# =========================
def fetch_organism_from_assembly(accession, max_retries=3):
    """
    单个 accession：nuccore → assembly → 取 Organism
    """
    for attempt in range(max_retries):
        try:
            # 1. 找对应的 Assembly
            handle = Entrez.elink(
                dbfrom="nuccore",
                db="assembly",
                id=accession
            )
            links = Entrez.read(handle)
            handle.close()
            
            if not links or not links[0]['LinkSetDb']:
                return None, None
            
            assembly_ids = [link['Id'] for link in links[0]['LinkSetDb'][0]['Link']]
            if not assembly_ids:
                return None, None
            
            # 2. 查 Assembly 摘要
            handle = Entrez.esummary(db="assembly", id=assembly_ids[0])
            asm_records = Entrez.read(handle)
            handle.close()
            
            asm = asm_records['DocumentSummarySet']['DocumentSummary'][0]
            organism = asm.get('Organism', '')
            taxid = asm.get('TaxId', None)
            
            return organism, taxid
        
        except Exception as e:
            print(f"⚠️ {accession} 第{attempt+1}次失败: {e}")
            time.sleep(2 + attempt * 2)
    return None, None

# ----------------------
# 执行：只对 unknown 重试
# ----------------------
print("\n🔍 处理 unknown 物种：accession → Assembly → ORGANISM")
unknown_mask = (df_ori['species_clean'] == 'unknown')
unknown_accs = df_ori.loc[unknown_mask, 'Accession_Number'].dropna().unique()
print(f"📌 待修复 unknown 数量：{len(unknown_accs)}")

fixed_count = 0

for idx, acc in enumerate(unknown_accs):
    if idx % 50 == 0:
        print(f"进度：{idx}/{len(unknown_accs)}")
    
    organism, taxid = fetch_organism_from_assembly(acc)
    
    if organism and organism.strip() != '':
        # 清洗
        clean_name, is_abn = clean_species(organism)
        
        # 更新全表
        mask = df_ori['Accession_Number'] == acc
        df_ori.loc[mask, 'species'] = organism
        df_ori.loc[mask, 'species_clean'] = clean_name
        df_ori.loc[mask, 'taxid'] = taxid
        df_ori.loc[mask, 'is_species_abnormal'] = is_abn
        fixed_count += 1
    
    time.sleep(0.5)

print(f"✅ 通过 Assembly 修复成功：{fixed_count}/{len(unknown_accs)}")

# 最终仍异常标记
df_ori['need_manual_check'] = (
    (df_ori['species_clean'] == 'unknown') | 
    (df_ori['is_species_abnormal'] == True)
)

# =========================
# 🧬 5. 构建 ID（不变）
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
# 🔗 6. 配对（不变）
# =========================
df_ori_small = df_ori[
    ['Accession_Number', 'ori_id', 'ori_start', 'ori_end', 'Intergenic_Sequence', 
     'species_clean', 'taxid', 'need_manual_check']
]
df_rip_small = df_rip[
    ['Accession_Number', 'rep_id', 'rep_start', 'rep_end', 'gene_id', 'product', 'translation']
]

print("🔗 正在进行 Ori–Rep 配对...")
df_merged = pd.merge(df_ori_small, df_rip_small, on='Accession_Number', how='left')
print(f"配对后行数: {len(df_merged)}")

# =========================
# 📏 7. 距离（不变）
# =========================
ori_mid = (df_merged['ori_start'] + df_merged['ori_end']) / 2
rep_mid = (df_merged['rep_start'] + df_merged['rep_end']) / 2
df_merged['distance'] = abs(ori_mid - rep_mid)

# =========================
# 📦 8. 输出（不变）
# =========================
df_final = pd.DataFrame({
    'OriC sequence': df_merged['Intergenic_Sequence'],
    'ori_id': df_merged['ori_id'],
    'plasmid_id': df_merged['Accession_Number'],
    
    'species': df_merged['species_clean'],
    'taxid': df_merged['taxid'],
    'need_manual_check': df_merged['need_manual_check'],
    
    'source': 'PLSDB',
    'pfamid_fast': np.nan,
    
    'rep_id': df_merged['rep_id'],
    'Rep_type_fast': df_merged['product'],
    'rep_seq': df_merged['translation'],
    
    'rep_dna_seq': np.nan,
    'full_replicon_seq': np.nan,
    'split': np.nan,
    
    'ori_start': df_merged['ori_start'],
    'ori_end': df_merged['ori_end'],
    'rep_start': df_merged['rep_start'],
    'rep_end': df_merged['rep_end'],
    'distance': df_merged['distance']
})

output_file = 'merged_final_optimized00000.csv'
df_final.to_csv(output_file, index=False)

print("\n✅ 全部完成！")
print(f"最终数据量: {len(df_final)}")
print(f"输出文件: {output_file}")

# =========================
# 📊 统计
# =========================
print("\n" + "="*60)
print("📊 最终物种质量校验")
print("="*60)
print(f"✅ 正常物种：{(df_final['species'] != 'unknown').sum()}")
print(f"⚠️  仍 unknown：{df_final['species'].eq('unknown').sum()}")
print(f"🧪 需人工检查：{df_final['need_manual_check'].sum()}")
print(f"✅ 物种正常率：{1 - df_final['species'].eq('unknown').mean():.2%}")
print("="*60)