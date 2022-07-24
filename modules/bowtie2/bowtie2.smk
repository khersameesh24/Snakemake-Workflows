rule run_alignment:
    """
    Align Reads to the Reference Genome
    Runs: Bowtie2
    """
    input:
        fastq_R1 = f"{out_dir}/Trimmed/{{sample}}_R1_001.fastq.gz",
        fastq_R2 = f"{out_dir}/Trimmed/{{sample}}_R2_001.fastq.gz",
    output:
        sam_out=f"{out_dir}/Aligned/{{sample}}",
    params:
        ref_genome=f"{in_dir}/Reference/hg19/hg19_karyotypicOrder_btw2_index",
    log:
        f"{out_dir}/logs/{{sample}}.bowtie2.log",
    benchmark:
        f"{out_dir}/benchmarks/{{sample}}.bowtie2.benchmarks.txt"
    threads: 8
    resources:
        mem_gb=8,
    shell:
        """bowtie2 -x {params.ref_genome} -1 {input.fastq_R1} -2 {input.fastq_R2} -S {output.sam_out} -p {threads}"""