import pandas as pd
import json

# ====================== 1. 筛选物种：只保留样本数 >5 的物种 ======================
df_species = pd.read_excel("csv_species_statistics.xlsx", sheet_name="Species_种统计")
valid_species = df_species[df_species["Count"] > 5]["Name"].tolist()
valid_species = sorted(list(set(valid_species)))

# ====================== 2. 构建词表 ======================
special_tokens = ["<PAD>", "<UNK_SPECIES>", "<SEP>", "<STOP>"]
aa_tokens = list("ACDEFGHIKLMNPQRSTVWY")
nt_tokens = list("atcg")

# 全部词汇
vocab = special_tokens + aa_tokens + nt_tokens + valid_species

# 映射
token2id = {tok: i for i, tok in enumerate(vocab)}
id2token = {i: tok for tok, i in token2id.items()}

# ====================== 3. 保存为词表文件 ======================
tokenizer_config = {
    "token2id": token2id,
    "id2token": id2token,
    "valid_species": valid_species,
    "vocab_size": len(vocab)
}

with open("origin_tokenizer.json", "w", encoding="utf-8") as f:
    json.dump(tokenizer_config, f, indent=2, ensure_ascii=False)

# ====================== 4. 最终信息 ======================
print(f"✅ 构建完成！")
print(f"总词汇量: {len(vocab)}")
print(f"有效物种数(>5): {len(valid_species)}")
print(f"词表已保存为: origin_tokenizer.json")