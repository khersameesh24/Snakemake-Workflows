
rule samtools_view:
    input:
        sam = "aligned/{sample}/Aligned.sam"
    output:
        bam = "aligned/{sample}/Aligned.bam"
    params: "--quiet"
    log:
        "logs/samtools/{sample}.log"
    threads: 8
    shell:
        """
        samtools view \
        -bS {input.sam} \
        --threads {threads} \
        -o {output.bam}
        """