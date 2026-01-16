#Load_environment
source ~/.bashrc
conda activate blast_env

#run_the_script(python name.py)
import os
import pandas as pd

#constants

input_folder = "/home/g89x126/Labwork_Bioinformatics/vfdb_csv_blast"

output_excel = "merged_vfdb_blast_results.xlsx"


#Collect_all_filtered_CSV_filenames

csv_files = [f for f in os.listdir(input_folder) if f.endswith("_filtered_results.csv")]

csv_files.sort()  # Sort filenames alphabetically



#Write_each_file_to_separate_sheet

with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:

    for filename in csv_files:

        try:

            file_path = os.path.join(input_folder, filename)

            df = pd.read_csv(file_path)



            #Extract_gene_name_from_sseqid

            df['gene'] = df['sseqid'].apply(lambda x: x.split('.')[0] if pd.notnull(x) else '')



            #Sheet_name=filename_before_'_filtered'

            sheet_name = filename.replace("_filtered_results.csv", "")

            if len(sheet_name) > 31:

                sheet_name = sheet_name[:31]


            df.to_excel(writer, sheet_name=sheet_name, index=False)

            print(f"Written sheet: {sheet_name}")

        except Exception as e:

            print(f"Error processing {filename}: {e}")

