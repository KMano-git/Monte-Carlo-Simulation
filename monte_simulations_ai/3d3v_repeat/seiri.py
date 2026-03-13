import glob
import re
import csv

# logファイルのリストを取得
log_files = glob.glob("result_10/run_*/run.log")

if not log_files:
    print("エラー: run.log ファイルが見つかりません。")
    exit()

# フォルダ名に含まれる数字（Run ID）でソートする関数
def get_run_id(filepath):
    # 例: "result/run_15/run.log" から 15 を抽出
    match = re.search(r'run_(\d+)', filepath)
    return int(match.group(1)) if match else 0

# ファイルリストをRun IDの昇順（1, 2, 3...30）に並び替え
log_files.sort(key=get_run_id)

# 出力用データのリスト（1行目にヘッダーを追加）
csv_data = [["Run_ID", "CL_CX", "CL_EL", "TL_CX", "TL_EL"]]

print(f"{len(log_files)} 個のログファイルを処理してCSVに出力します...")

for file_path in log_files:
    run_id = get_run_id(file_path)
    cl_cx = cl_el = tl_cx = tl_el = None
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(" CL:"):
                parts = line.split()
                # parts[0]="CL:", parts[1]=EI, parts[2]=CX, parts[3]=EL
                cl_cx = float(parts[2])
                cl_el = float(parts[3])
            elif line.startswith(" TL:"):
                parts = line.split()
                # parts[0]="TL:", parts[1]=EI, parts[2]=CX, parts[3]=EL
                tl_cx = float(parts[2])
                tl_el = float(parts[3])
    
    # 取得したデータをリストに追加
    csv_data.append([run_id, cl_cx, cl_el, tl_cx, tl_el])

# CSVファイルへの書き込み
output_filename = "ensemble_results_10.csv"
with open(output_filename, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerows(csv_data)

print(f"完了しました！抽出した生データを '{output_filename}' に保存しました。")