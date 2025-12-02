from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import os
import pandas as pd
import logging

app = Flask(__name__)
CORS(app)

logging.basicConfig(level=logging.DEBUG)

# 甲基化文件目录
METH_BASE_DIR = r"D:\4four_year\20241203encode\chr_1k_region_csv"

# 自动扫描目录，生成所有 ChrXX[AB]-Context 的映射
methylation_mapping = {}

for fname in os.listdir(METH_BASE_DIR):
    # 只处理甲基化结果文件
    if fname.endswith(".methylationlevel.final.csv"):
        # 例如 fname = "Chr01A_CG.methylationlevel.final.csv"
        # 先取去掉后缀的前半部分：Chr01A_CG
        prefix = fname.split(".")[0]          # "Chr01A_CG"
        try:
            chr_name, context = prefix.split("_")  # "Chr01A", "CG/CHG/CHH"
        except ValueError:
            # 以防万一遇到奇怪文件名，跳过
            continue

        key = f"{chr_name}-{context}"         # 比如 "Chr01A-CG"
        full_path = os.path.join(METH_BASE_DIR, fname)
        methylation_mapping[key] = full_path

# 看看加载了哪些 key（可选调试）
print("Loaded methylation files:")
for k, v in methylation_mapping.items():
    print("  ", k, "->", v)

GENE_ANNOTATION_FILE = r"D:\4four_year\result\HM\result\peak_csv_modified\Cenchrus_fungigraminus_new.csv"  # <-- 这里你要填真实路径

@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    response.headers.add('X-Content-Type-Options', 'nosniff')
    return response


@app.route("/")
def home():
    return send_from_directory(os.path.dirname(os.path.abspath(__file__)), "methylation.html")


# ------------------------------
# 1️⃣ 区间搜索
# ------------------------------
@app.route("/search_methylation")
def search_methylation():
    try:
        chr_query = request.args.get("chr")
        start = request.args.get("start")
        end = request.args.get("end")
        context = request.args.get("context", "CG")

        app.logger.info(f"methylation query: chr={chr_query}, start={start}, end={end}, context={context}")

        if not chr_query or not start or not end:
            return jsonify({"error": "Missing required parameters (chr/start/end)"}), 400

        file_key = f"{chr_query}-{context}"
        file_path = methylation_mapping.get(file_key)

        if not file_path:
            return jsonify({"error": f"No methylation file configured for {file_key}"}), 404

        df = pd.read_csv(file_path, keep_default_na=False)
        start = int(start)
        end = int(end)

        filtered = df[
            (df["chr"] == chr_query) &
            (df["start"] <= end) &
            (df["end"] >= start)
        ].copy()

        filtered = filtered.sort_values(by="start")
        filtered = filtered.where(pd.notnull(filtered), None)

        return jsonify(filtered.to_dict(orient="records"))

    except Exception as e:
        app.logger.error(str(e), exc_info=True)
        return jsonify({"error": f"Server error: {str(e)}"}), 500


# ------------------------------
# 2️⃣ 基因搜索
# ------------------------------
@app.route("/search_gene_methylation")
def search_gene_methylation():
    try:
        gene_id = request.args.get("gene_id")
        context = request.args.get("context", "CG")

        if not gene_id:
            return jsonify({"error": "Missing gene_id parameter"}), 400

        if not os.path.exists(GENE_ANNOTATION_FILE):
            return jsonify({"error": "Gene annotation file not found"}), 500

        gene_df = pd.read_csv(GENE_ANNOTATION_FILE, keep_default_na=False)
        gene_row = gene_df[gene_df["Gene_ID"].astype(str).str.strip() == gene_id.strip()]

        if gene_row.empty:
            return jsonify({"error": f"Gene ID '{gene_id}' not found"}), 404

        gene = gene_row.iloc[0]
        gene_chr = str(gene["Chromosome"])
        gene_start = int(gene["Start"])
        gene_end = int(gene["End"])

        file_key = f"{gene_chr}-{context}"
        file_path = methylation_mapping.get(file_key)

        if not file_path:
            return jsonify({"error": f"No methylation file configured for {file_key}"}), 404

        df = pd.read_csv(file_path, keep_default_na=False)

        mask = (
            (df["chr"] == gene_chr) &
            (df["start"] <= gene_end) &
            (df["end"] >= gene_start)
        )

        filtered = df[mask].copy().sort_values(by="start")
        filtered = filtered.where(pd.notnull(filtered), None)

        return jsonify({
            "gene_info": {
                "gene_id": gene_id,
                "chr": gene_chr,
                "start": gene_start,
                "end": gene_end,
                "description": gene.get("Description", "")
            },
            "windows": filtered.to_dict(orient="records")
        })

    except Exception as e:
        app.logger.error(str(e), exc_info=True)
        return jsonify({"error": f"Server error: {str(e)}"}), 500


if __name__ == "__main__":
    app.run(debug=True, port=5001)
