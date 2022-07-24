rule run_bedtools:
    """
    Get Coverage from the bam file
    indexed in previous step
    Runs: Bedtools
    """
    input:
        sorted_bam=f"{out_dir}/Aligned/{{sample}}.sorted.bam",
    output:
        bed=f"{out_dir}/Bedfiles/{{sample}}.bed",
    params:
        ref_bed=f"{in_dir}/Reference/hg19/hg19_refseq.bed",
    log:
        f"{out_dir}/logs/{{sample}}.bedtools.log",
    benchmark:
        f"{out_dir}/benchmarks/{{sample}}.bedtools.benchmarks.txt"
    threads: 8
    resources:
        mem_gb=8,
    shell:
        """bedtools multicov -bams {input.sorted_bam} -bed {params.ref_bed} > {output.bed}"""