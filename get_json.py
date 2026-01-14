import pandas as pd
import json
import os
import re

# 蛋白质序列（A链）
PROTEIN_SEQUENCE = "IEKKKSFAKGMGVKSTLVSGSKVYMTTFAEGSDARLEKIVEGDSIRSVNEGEAFSAEMADKNAGYKIGNAKFSHPKGYAVVANNPLYTGPVQQDMLGLKETLEKRYFGESADGNDNICIQVIHNILDIEKILAEYITNAAYAVNNISGLDKDIIGFGKFSTVYTYDEFKDPEHHRAAFNNNDKLINAIKAQYDEFDNFLDNPRLGYFGQAFFSKEGRNYIINYGNECYDILALLSGLAHWVVANNEEESRISRTWLYNLDKNLDNEYISTLNYLYDRITNELTNSFSKNSAANVNYIAETLGINPAEFAEQYFRFSIMKEQKNLGFNITKLREVMLDRKDMSEIRKNHKVFDSIRTKVYTMMDFVIYRYYIEEDAKVAAANKSLPDNEKSLSEKDIFVINLRGSFNDDQKDALYYDEANRIWRKLENIMHNIKEFRGNKTREYKKKDAPRLPRILPAGRDVSAFSKLMYALTMFLDGKEINDLLTTLINKFDNIQSFLKVMPLIGVNAKFVEEYAFFKDSAKIADELRLIKSFARMGEPIADARRAMYIDAIRILGTNLSYDELKALADTFSLDENGNKLKKGKHGMRNFIINNVISNKRFHYLIRYGDPAHLHEIAKNEAVVKFVLGRIADIQKKQGQNGKNQIDRYYETCIGKDKGKSVSEKVDALTKIITGMNYDQFDKKRSVIEDTGRENAEREKFKKIISLYLTVIYHILKNIVNINARYVIGFHCVERDAQLYKEKGYDINLKKLEEKGFSSVTKLCAGIDETAPDKRKDVEKEMAERAKESIDSLESANPKLYANYIKYSDEKKAEEFTRQINREKAKTALNAYLRNTKWNVIIREDLLRIDNKTCTLFANKAVALEVARYVHAYINDIAEVNSYFQLYHYIMQRIIMNERYEKSSGKVSEYFDAVNDEKKYNDRLLKLLCVPFGYCIPRFKNLSIEALFDRNEAAKFDKEKKKVSGNS"

# 修正后的输入文件列表（DM替换第二个PM）
input_files = [
    "guide_sequences_CD.fasta_filtered_by_seq.csv",
    "guide_sequences_CI.fasta_filtered_by_seq.csv",
    "guide_sequences_DD.fasta_filtered_by_seq.csv",
    "guide_sequences_DI.fasta_filtered_by_seq.csv",
    "guide_sequences_PM.fasta_filtered_by_seq.csv",
    "guide_sequences_DM.fasta_filtered_by_seq.csv",
    "guide_sequences_RDM.fasta_filtered_by_seq.csv",
    "guide_sequences_RTM.fasta_filtered_by_seq.csv",
    "guide_sequences_SD.fasta_filtered_by_seq.csv",
    "guide_sequences_SI.fasta_filtered_by_seq.csv",
    "guide_sequences_SM.fasta_filtered_by_seq.csv",
    "guide_sequences_TM.fasta_filtered_by_seq.csv"
]

def clean_guide_id(guide_id):
    """清理guide_id以作为文件名"""
    if pd.isna(guide_id):
        return ""
    
    # 转换为字符串
    guide_id_str = str(guide_id)
    
    # 移除非法字符，保留字母、数字、下划线、连字符、点、冒号和管道符
    # 先替换冒号和管道符为下划线
    cleaned = guide_id_str.replace(':', '_').replace('|', '_')
    
    # 移除其他非法字符
    cleaned = re.sub(r'[<>:"/\\|?*]', '_', cleaned)
    
    # 确保不以点或空格开头/结尾
    cleaned = cleaned.strip('. ')
    
    return cleaned

def read_csv_with_auto_delimiter(file_path):
    """尝试不同的分隔符读取CSV文件"""
    # 先尝试逗号分隔
    try:
        df = pd.read_csv(file_path)
        if len(df.columns) > 1:
            print(f"  使用逗号分隔符，读取到 {len(df)} 行，{len(df.columns)} 列")
            return df
    except:
        pass
    
    # 尝试制表符分隔
    try:
        df = pd.read_csv(file_path, sep='\t')
        if len(df.columns) > 1:
            print(f"  使用制表符分隔符，读取到 {len(df)} 行，{len(df.columns)} 列")
            return df
    except:
        pass
    
    # 尝试分号分隔
    try:
        df = pd.read_csv(file_path, sep=';')
        if len(df.columns) > 1:
            print(f"  使用分号分隔符，读取到 {len(df)} 行，{len(df.columns)} 列")
            return df
    except:
        pass
    
    # 如果都不行，尝试读取第一行判断
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline()
        
        if ',' in first_line:
            df = pd.read_csv(file_path)
            print(f"  检测到逗号，使用逗号分隔符")
        elif '\t' in first_line:
            df = pd.read_csv(file_path, sep='\t')
            print(f"  检测到制表符，使用制表符分隔符")
        else:
            # 默认使用逗号
            df = pd.read_csv(file_path)
            print(f"  使用默认逗号分隔符")
        
        return df
    except Exception as e:
        raise Exception(f"无法读取文件 {file_path}: {str(e)}")

def process_csv_file(csv_file_path, output_dir):
    """
    处理单个CSV文件，提取序列并生成JSON文件
    """
    try:
        print(f"\n正在处理文件: {csv_file_path}")
        
        # 检查文件是否存在
        if not os.path.exists(csv_file_path):
            print(f"错误: 文件不存在: {csv_file_path}")
            return 0
        
        # 读取CSV文件（自动检测分隔符）
        df = read_csv_with_auto_delimiter(csv_file_path)
        
        print(f"成功读取 {len(df)} 行数据")
        
        # 检查必要的列是否存在
        required_columns = ['full_target_sequence', 'guide_seq', 'guide_id']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"错误: 文件 {csv_file_path} 缺少以下列: {missing_columns}")
            print(f"可用的列: {list(df.columns)}")
            return 0
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        json_count = 0
        skipped_count = 0
        
        # 处理每一行
        for index, row in df.iterrows():
            # 获取序列和guide_id
            target_sequence = row['full_target_sequence']
            guide_sequence = row['guide_seq']
            guide_id = row['guide_id']
            
            # 检查序列是否有效
            if pd.isna(target_sequence) or pd.isna(guide_sequence):
                print(f"警告: 行 {index+1} 的序列为空，跳过")
                skipped_count += 1
                continue
            
            # 清理guide_id作为文件名
            safe_guide_id = clean_guide_id(guide_id)
            if not safe_guide_id:
                safe_guide_id = f"sequence_{index+1}"
            
            # 创建JSON结构
            json_data = {
                "name": safe_guide_id,
                "sequences": [
                    {
                        "protein": {
                            "id": ["A"],
                            "sequence": PROTEIN_SEQUENCE
                        }
                    },
                    {
                        "rna": {
                            "id": ["B"],
                            "sequence": str(target_sequence)
                        }
                    },
                    {
                        "rna": {
                            "id": ["C"],
                            "sequence": str(guide_sequence)
                        }
                    }
                ],
                "modelSeeds": [1],
                "dialect": "alphafold3",
                "version": 1
            }
            
            # 保存JSON文件
            json_filename = f"{safe_guide_id}.json"
            json_path = os.path.join(output_dir, json_filename)
            
            with open(json_path, 'w', encoding='utf-8') as json_file:
                json.dump(json_data, json_file, indent=2, ensure_ascii=False)
            
            json_count += 1
            
            # 显示进度
            if (index + 1) % 100 == 0:
                print(f"  已处理 {index+1}/{len(df)} 行...")
        
        print(f"处理完成: {csv_file_path}")
        print(f"  -> 输出目录: {output_dir}")
        print(f"  -> 生成JSON文件: {json_count} 个")
        if skipped_count > 0:
            print(f"  -> 跳过无效行: {skipped_count} 个")
        
        return json_count
        
    except Exception as e:
        print(f"处理文件 {csv_file_path} 时出错: {str(e)}")
        import traceback
        traceback.print_exc()
        return 0

def main():
    print("开始处理CSV文件并生成AlphaFold3 JSON文件...")
    print("=" * 60)
    
    total_json_count = 0
    
    # 为每个CSV文件创建输出目录并处理
    for csv_file in input_files:
        # 生成输出目录名称（基于CSV文件名）
        base_name = os.path.splitext(csv_file)[0]
        output_dir = f"json_output_{base_name}"
        
        # 处理文件
        json_count = process_csv_file(csv_file, output_dir)
        total_json_count += json_count
        
        print("-" * 50)
    
    print("=" * 60)
    print(f"所有文件处理完成!")
    print(f"总共生成 {total_json_count} 个JSON文件")
    
    # 显示输出目录结构
    print("\n输出目录结构:")
    for csv_file in input_files:
        base_name = os.path.splitext(csv_file)[0]
        output_dir = f"json_output_{base_name}"
        if os.path.exists(output_dir):
            file_count = len([f for f in os.listdir(output_dir) if f.endswith('.json')])
            print(f"  {output_dir}/ - {file_count} 个JSON文件")

# 单个文件处理函数（如果需要单独处理某个文件）
def process_single_file(csv_file_path):
    """处理单个CSV文件"""
    if not os.path.exists(csv_file_path):
        print(f"错误: 文件不存在: {csv_file_path}")
        return 0
    
    # 从文件名生成输出目录
    base_name = os.path.splitext(os.path.basename(csv_file_path))[0]
    output_dir = f"json_output_{base_name}"
    
    json_count = process_csv_file(csv_file_path, output_dir)
    return json_count

if __name__ == "__main__":
    # 运行主函数处理所有文件
    main()
    
    # 如果只想处理单个文件，可以使用下面的代码：
    # json_count = process_single_file("guide_sequences_CD.fasta_filtered_by_seq.csv")
    # print(f"生成了 {json_count} 个JSON文件")