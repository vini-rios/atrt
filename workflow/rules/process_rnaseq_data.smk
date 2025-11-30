import os
import glob
import pathlib
import zipfile


configfile: "config.yaml"


include: "common.smk"


PROJECT = ["PRJNA565377"]
OUTPUT_DIR = "data/"


rule get_fastq_list:
    output:
        os.path.join(OUTPUT_DIR, "{project}.files_to_download.tsv"),
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),  # Infers directory from the output file
        project=lambda wildcards: wildcards.project,
    log:
        os.path.join("logs", "get_fastq_list_{project}.log"),
    conda:
        "../envs/core.yml"
    shell:
        """
        scripts/get_fastq_ena.sh -l {params.project} --outputdir {params.output_dir}
        """


checkpoint download_fastq:
    input:
        os.path.join(OUTPUT_DIR, "{project}.files_to_download.tsv"),
    output:
        os.path.join(OUTPUT_DIR, "raw", "{project}_download_complete.flag"),
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),  # Infers directory from the output file
        project=lambda wildcards: wildcards.project,
    log:
        os.path.join("logs", "download_fastq_{project}.log"),
    conda:
        "../envs/core.yml"
    shell:
        """
        scripts/get_fastq_ena.sh --retry {params.project} --outputdir {params.output_dir}
        touch {output}
        """


rule check_download:
    input:
        flag=os.path.join(OUTPUT_DIR, "raw", "{project}_download_complete.flag"),
        file_list=os.path.join(OUTPUT_DIR, "{project}.files_to_download.tsv"),
    output:
        os.path.join(OUTPUT_DIR, "raw", "{project}_all_files_verified.flag"),
    log:
        os.path.join("logs", "check_download_{project}.log"),
    conda:
        "../envs/core.yml"
    script:
        "../scripts/check_download.py"


rule fastp_paired_end:
    input:
        flag=os.path.join(OUTPUT_DIR, "raw", "{project}_all_files_verified.flag"),
        r1="data/raw/{project}__{base}_1.fastq.gz",
        r2="data/raw/{project}__{base}_2.fastq.gz",
    output:
        r1="results/fastp/{project}/{base}_1_trimmed.fastq.gz",
        r2="results/fastp/{project}/{base}_2_trimmed.fastq.gz",
        html="results/fastp/{project}/{base}_fastp.html",
        json="results/fastp/{project}/{base}_fastp.json",
    log:
        os.path.join("logs", "{project}", "fastp_{base}.log"),
    conda:
        "../envs/fastp.yml"
    threads: 4
    shell:
        """
        mkdir -p results/fastp/{wildcards.project}
        fastp -i {input.r1} \
              -I {input.r2} \
              -o {output.r1} \
              -O {output.r2} \
              --html {output.html} \
              --json {output.json} \
              --thread {threads} \
              --detect_adapter_for_pe
        """


rule fastp_single_end:
    input:
        flag=os.path.join(OUTPUT_DIR, "raw", "{project}_all_files_verified.flag"),
        r1="data/raw/{project}__{base}.fastq.gz",
    output:
        r1="results/fastp/{project}/{base}_trimmed.fastq.gz",
        html="results/fastp/{project}/{base}_fastp.html",
        json="results/fastp/{project}/{base}_fastp.json",
    log:
        os.path.join("logs", "{project}", "fastp_{base}.log"),
    conda:
        "../envs/fastp.yml"
    threads: 4
    shell:
        """
        mkdir -p results/fastp/{wildcards.project}
        fastp -i {input.r1} \
              -o {output.r1} \
              --html {output.html} \
              --json {output.json} \
              --thread {threads}
        """


ruleorder: fastp_paired_end > fastp_single_end


rule fastp_project_complete:
    input:
        get_all_fastp_outputs,
    output:
        "results/fastp/{project}_fastp_complete.flag",
    log:
        os.path.join("logs", "{project}", "fastp_project_complete.log"),
    conda:
        "../envs/core.yml"
    shell:
        """
        touch {output}
        """


rule download_references:
    output:
        transcriptome="references/homo_sapiens.transcriptome.fa.gz",
        genome="references/homo_sapiens.genome.fa.gz",  # For decoys
        gtf="references/homo_sapiens.gtf.gz",
    params:
        ref_ver=config["reference_version"],
    log:
        os.path.join("logs", "download_reference.log"),
    conda:
        "../envs/core.yml"
    shell:
        """
        mkdir -p references 
        # Download from Ensembl
        wget --timestamping -O {output.transcriptome} https://ftp.ensembl.org/pub/release-{params.ref_ver}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        wget --timestamping -O {output.genome} https://ftp.ensembl.org/pub/release-{params.ref_ver}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        wget --timestamping -O {output.gtf} https://ftp.ensembl.org/pub/release-{params.ref_ver}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{params.ref_ver}.gtf.gz
        """


rule build_salmon_index:
    input:
        transcriptome="references/homo_sapiens.transcriptome.fa.gz",
        genome="references/homo_sapiens.genome.fa.gz",
    output:
        index=directory("references/salmon_index"),
    log:
        os.path.join("logs", "build_salmon_index.log"),
    conda:
        "../envs/salmon.yml"
    threads: 4
    shell:
        """
        # Extract decoy names from genome
        grep "^>" <(gunzip -c {input.genome}) | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt
        
        # Combine transcriptome and genome for decoy-aware index
        cat {input.transcriptome} {input.genome} > transcriptome_genome.fa.gz
        
        salmon index -t transcriptome_genome.fa.gz -d decoys.txt -i {output.index} --threads {threads}
        
        # Cleanup
        rm decoys.txt transcriptome_genome.fa.gz
        """


rule salmon_quant_single:
    input:
        flag="results/fastp/{project}_fastp_complete.flag",
        r1="results/fastp/{project}/{base}_trimmed.fastq.gz",
        index="references/salmon_index",
    output:
        "results/quantification/salmon/{project}/{base}/quant.sf",
    log:
        os.path.join("logs", "{project}/{base}" "salmon_quant_single.log"),
    conda:
        "../envs/salmon.yml"
    threads: 4
    shell:
        """
        salmon quant --gcBias -i {input.index} -l A -r {input.r1} \
            -o results/quantification/salmon/{wildcards.project}/{wildcards.base} --threads {threads}
        """


rule salmon_quant_paired:
    input:
        flag="results/fastp/{project}_fastp_complete.flag",
        r1="results/fastp/{project}/{base}_1_trimmed.fastq.gz",
        r2="results/fastp/{project}/{base}_2_trimmed.fastq.gz",
        index="references/salmon_index",
    output:
        "results/quantification/salmon/{project}/{base}/quant.sf",
    log:
        os.path.join("logs", "{project}/{base}" "salmon_quant_paired.log"),
    conda:
        "../envs/salmon.yml"
    threads: 4
    shell:
        """
        salmon quant --gcBias -i {input.index} -l A -1 {input.r1} -2 {input.r2} \
            -o results/quantification/salmon/{wildcards.project}/{wildcards.base} --threads {threads}
        """


ruleorder: salmon_quant_paired > salmon_quant_single


rule get_metadata_from_ena:
    input:
        os.path.join(OUTPUT_DIR, "{project}.files_to_download.tsv"),
    output:
        os.path.join(OUTPUT_DIR, "{project}.metadata"),
    conda:
        "../envs/core.yml"
    shell:
        """
        wget --output-document="{output}" "https://www.ebi.ac.uk/ena/portal/api/links/study?accession={wildcards.project}&result=read_run"
        """


rule clean_metadata:
    input:
        os.path.join(OUTPUT_DIR, "{project}.metadata"),
    output:
        os.path.join(OUTPUT_DIR, "{project}.metadata.clean"),

    shell:
        """
        ./scripts/format_metadata.r {wildcards.project}
        """
