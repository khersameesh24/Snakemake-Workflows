rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("ref/star_genome"),
    threads: 8
    params:
        gtf_file = "--sjdbGTFfile resources/genome.gtf",
        overhang = "--sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    shell:
        """
        STAR \
        --runMode GenomeGenerate \
        --genomeDir ref/ \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFFile {input.annotation} \
        --runThreadN {threads}
        """
