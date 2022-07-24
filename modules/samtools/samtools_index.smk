rule samtools_index:
    input:
        bam = "aligned/{sample}/Aligned.bam"
    output:
        sorted_bam = "aligned/{sample}/Aligned.Sorted.bam"
    params: "--quiet"
    log:
        "logs/samtools/{sample}.log"
    threads: 8
    shell:
        """
        samtools sort \
        {input.bam} \
        -o {output.sorted_bam} \
        -@ {threads}
        """