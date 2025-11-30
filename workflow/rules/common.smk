import os
import glob


def get_samples_and_types(project):
    raw_dir = os.path.join(OUTPUT_DIR, "raw")
    samples = {}
    # Look for all fastq.gz files for this project
    pattern = os.path.join(raw_dir, f"{project}__*.fastq.gz")
    all_files = glob.glob(pattern)
    for file_path in all_files:
        filename = os.path.basename(file_path)
        base_name = filename.replace(f"{project}__", "").replace(".fastq.gz", "")
        if base_name.endswith("_1"):
            sample_base = base_name[:-2]
            if sample_base not in samples:
                samples[sample_base] = {"type": None, "r1": None, "r2": None}
            samples[sample_base]["r1"] = file_path
            samples[sample_base]["type"] = "paired"
        elif base_name.endswith("_2"):
            sample_base = base_name[:-2]
            if sample_base not in samples:
                samples[sample_base] = {"type": None, "r1": None, "r2": None}
            samples[sample_base]["r2"] = file_path
            samples[sample_base]["type"] = "paired"
        else:
            samples[base_name] = {
                "type": "single",
                "r1": file_path,
                "r2": None
            }
    return samples

def get_all_fastp_outputs(wildcards):
    samples = get_samples_and_types(wildcards.project)
    outputs = []
    for sample, info in samples.items():
        if info["type"] == "paired":
            outputs.extend([
                f"results/fastp/{wildcards.project}/{sample}_1_trimmed.fastq.gz",
                f"results/fastp/{wildcards.project}/{sample}_2_trimmed.fastq.gz",
                #f"results/fastp/{wildcards.project}/{sample}_fastp.html",
                #f"results/fastp/{wildcards.project}/{sample}_fastp.json"
            ])
        else:
            outputs.extend([
                f"results/fastp/{wildcards.project}/{sample}_trimmed.fastq.gz",
                #f"results/fastp/{wildcards.project}/{sample}_fastp.html",
                #f"results/fastp/{wildcards.project}/{sample}_fastp.json"
            ])
    return outputs



def get_all_quantification_outputs(project):
    samples = get_samples_and_types(project)
    outputs = []
    
    for sample, info in samples.items():
        outputs.append(f"results/quantification/salmon/{project}/{sample}/quant.sf")
    return outputs

