rule fastqc:
    input:
        input_path = "reads/{sample}.fastq"
    output:
        html="qc/fastqc/{sample}.html",
        output_path = "qc/fastqc/{sample}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 8
    shell:
        """
        fastqc \
        {snakemake.params} \
        -t {snakemake.threads} \
        --outdir {output.output_path} \
        {input.input_path}"
        " {log}"
        """