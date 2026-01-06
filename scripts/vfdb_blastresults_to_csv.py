
#Load_environment
source ~/.bashrc
conda activate blast_env

import os
import pandas as pd

#Constants
MIN_PIDENT = 50  #Minimum_percent_identity_threshold

OUTPUT_FOLDER = "/home/username/Labwork_Bioinformatics/vfdb_csv_blast"  #Output_folder_for_CSV_files

MAIN_FOLDER = "/home/username/Labwork_Bioinformatics/BLAST_results_vfdb"  #The_main_folder_containing_strain_folders

#subfolder_in_the_main_folder

for strain_folder in os.listdir(MAIN_FOLDER):

    folder_path = os.path.join(MAIN_FOLDER, strain_folder)


    #Only_process_subfolders
    if os.path.isdir(folder_path):

        #filtered_results_text_file

        filtered_file = os.path.join(folder_path, f"{strain_folder}_filtered_results.txt")

        

        if os.path.isfile(filtered_file):

            try:

                #Read_BLAST_results_into_a_DataFrame

                df = pd.read_csv(filtered_file, sep="\t", header=None)

                df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",

                              "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
              

                #Convert_pident_column_to_numeric

                df['pident'] = pd.to_numeric(df['pident'], errors='coerce')
              

                #best_hit_per_query

                df.sort_values(["qseqid", "evalue"], ascending=[True, True], inplace=True)

                best_hits = df.drop_duplicates(subset="qseqid", keep="first")


                #Apply_identity_filter

                filtered_df = best_hits[best_hits["pident"] > MIN_PIDENT]
              

                #Define_output_CSV_file_path

                output_csv = os.path.join(OUTPUT_FOLDER, f"{strain_folder}_filtered_results.csv")
                

                #filtered_DataFrame_to_CSV_file

                filtered_df.to_csv(output_csv, index=False)


                print(f"Processed {strain_folder} and saved to {output_csv}")

            except Exception as e:

                print(f"Error processing {strain_folder}: {e}")

        else:

            print(f"File {filtered_file} does not exist or is not a valid file.")

