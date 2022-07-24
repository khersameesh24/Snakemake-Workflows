import pathlib

def get_bam_files(bam_path: Path):
    bam_files = [files for files in bam_files if files.endswith("Aligned.Sorted.bam") ]
    return bam_files

rule feature_counts:
    input:
        path = get_bam_files,
        annotation="resources/genome.gtf"
    output:
        matrix = "counts/counts.csv"
    params:
    log:
        "logs/featureCounts.log"
    threads: 8
    shell:
        """
        featureCounts \
        -a {input.annotation} \
        -o {output.matrix} \
        -T {threads} \
        {input.path}/*.bam
        """