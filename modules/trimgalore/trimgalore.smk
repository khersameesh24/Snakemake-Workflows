rule run_trimgalore:
    """
    Trim Adapter Sequences from 
    Sequenced Reads
    Runs: TrimGalore
    """
    input:
        fastq_R1=f"{in_dir}/{{sample}}_R1_001.fastq.gz",
        fastq_R2=f"{in_dir}/{{sample}}_R2_001.fastq.gz",
    params:
        outdir=f"{out_dir}/Trimmed",
    log:
        f"{out_dir}/logs/{{sample}}.trimgalore.log",
    benchmark:
        f"{out_dir}/benchmarks/{{sample}}.trimgalore.benchmarks.txt"
    threads: 8
    resources:
        mem_gb=8,
    shell:
        """trim_galore --fastqc --paired {input.fastq_R1} {input.fastq_R2} -o {params.outdir}"""
