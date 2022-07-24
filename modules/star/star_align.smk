rule star_pe_multi:
    input:
        fastq = ["reads/{sample}_R1.fastq", "reads/{sample}_R2.fastq"]
    output:
        sam = "aligned/{sample}/Aligned.out.sam",
        log = "aligned/{sample}/Log.out",
    log:
        "logs/star/pe/{sample}.log",
    params:
        index = "ref/star_genome",
        tmpdir = "tmpdir/temp"
    threads: 8
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input.fastq} \
        --outTmpDir {params.tmpdir} \
        --outFileNamePrefix {output.sam}/ \
        --outStd {output.log}
        """