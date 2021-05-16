run_folder = config["run_folder"]
fastq_ext = "001.fastq.gz"
in_dir = f"/home/JNCASR-VBL/{run_folder}"
out_dir = f"/home/JNCASR-VBL/{run_folder}/Analysis"

all_terminal_file = []


rule all:
    input:
        all_terminal_file,


# Run FastQC to check read quality
rule run_fastqc:
    """Check Quality of the sequenced 
        Reads
    Runs: FastQC
    """
    input:
        fastq_R1=f"{in_dir}/{{sample}}_R1_{fastq_ext}",
        fastq_R2=f"{in_dir}/{{sample}}_R2_{fastq_ext}",
    params:
        outdir=out_dir,
    log:
        f"{out_dir}/logs/{{sample}}_{fastq_ext}.fastqc.log",
    benchmark:
        f"{out_dir}/benchmarks/{{sample}}_{fastq_ext}.fastqc.benchmarks.txt"
    threads: 8
    resources:
        mem_gb=8,
    run:
        shell("""fastqc {input.fastq_R1} --outdir={params.outdir} -t {threads}""")
        shell("""fastqc {input.fastq_R2} --outdir={params.outdir} -t {threads}""")


# Run TrimGalore to trim adapter sequences
rule run_trimgalore:
    """Trim Adapter Sequences from 
        Sequenced Reads
    Runs: TrimGalore
    """
    input:
        fastq_R1=f"{in_dir}/{{sample}}_R1_{fastq_ext}",
        fastq_R2=f"{in_dir}/{{sample}}_R2_{fastq_ext}",
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


# Run alignment to align reads onto ref genome
rule run_alignment:
    """Align Reads to the Reference
        Genome
    Runs: Bowtie2
    """
    input:
        fastq_R1=f"{out_dir}/Trimmed/{{sample}}_R1_{fastq_ext}",
        fastq_R2=f"{out_dir}/Trimmed/{{sample}}_R2_{fastq_ext}",
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


# Run Samtools view to convert sam to bam
rule run_sam_to_bam:
    """Convert sam files from previous 
        step to bam format
    Runs: Samtools view
    """
    input:
        sam=f"{out_dir}/Aligned/{{sample}}",
    output:
        bam=f"{out_dir}/Aligned/{{sample}}.bam",
    log:
        f"{out_dir}/logs/{{sample}}.samtools_view.log",
    benchmark:
        f"{out_dir}/benchmarks/{{sample}}.samtools_view.benchmarks.txt"
    threads: 8
    resources:
        mem_gb=8,
    shell:
        """samtools view -bS {input.sam} -o {output.bam} -@ {threads}"""


# Run Samtools sort to sort bam files
rule run_sort_bam:
    """Sort bam files from 
        previous step
    Runs: Samtools sort
    """
    input:
        bam=f"{out_dir}/Aligned/{{sample}}.bam",
    output:
        sorted_bam=f"{out_dir}/Aligned/{{sample}}.sorted.bam",
    log:
        f"{out_dir}/logs/{{sample}}.samtools_sort.log",
    benchmark:
        f"{out_dir}/benchmarks/{{sample}}.samtools_sort.benchmarks.txt"
    threads: 8
    resources:
        mem_gb=8,
    shell:
        """samtools sort {input.bam} -o {output.sorted_bam}"""


# Run Samtools index to index the sorted bam file
rule run_index_bam:
    """Index the sorted bam files 
        from previous step
    Runs : Samtools index
    """
    input:
        sorted_bam=f"{out_dir}/Aligned/{{sample}}.sorted.bam",
    log:
        f"{out_dir}/logs/{{sample}}.samtools_index.log",
    benchmark:
        f"{out_dir}/benchmarks/{{sample}}.samtools_index.benchmarks.txt"
    threads: 8
    resources:
        mem_gb=8,
    shell:
        """samtools index {input.sorted_bam}"""


# Run Bedtools to get Coverage
rule run_bedtools:
    """Get Coverage from the bam file
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
