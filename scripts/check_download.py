import os
import pathlib

input_file_list = snakemake.input.file_list
output_flag = snakemake.output[0]
project = snakemake.wildcards.project
output_dir = os.path.dirname(output_flag)

with open(input_file_list) as f:
    lines = f.readlines()[:]
    expected_files = [f"{project}__" + str(os.path.basename(line.strip().split("\t")[1])) for line in lines]

missing_files = []
for filename in expected_files:
    filepath = os.path.join("data/raw/", filename) 
    if not os.path.exists(filepath):
        missing_files.append(filename)

if missing_files:
    raise ValueError(f"Missing files after download: {missing_files}")
else:
    print("All files downloaded successfully.")
    pathlib.Path.touch(output_flag)