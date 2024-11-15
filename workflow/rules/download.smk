

rule download_fastq_pe:
    output:
        temp("data/pe/{accession}_1.fastq.gz"),
        temp("data/pe/{accession}_2.fastq.gz"),
    log:
        "results/logs/download_fastq/pe/{accession}.log",
    params:
        out_prefix="data/pe",
        extra=config["download_fastq_pe"]["params"],
        tempdir=lambda wildcards: "data/pe/{}.fasterq".format(wildcards.accession),
    threads: config["download_fastq_pe"]["threads"]
    resources:
        max_downloads=1
    retries: config["download_fastq_pe"]["retries"]
    conda:
        "../envs/sra-tools.yaml"
    shell:
        "("
        "fasterq-dump"
        " --split-3 --skip-technical --force"
        " --threads {threads}"
        " --outdir {params.out_prefix}"
        " --temp {params.tempdir}"
        " {params.extra}"
        " {wildcards.accession}"
        "; pigz --processes {threads} {params.out_prefix}/{wildcards.accession}_1.fastq"
        "; pigz --processes {threads} {params.out_prefix}/{wildcards.accession}_2.fastq"
        "; rm -fr {params.tempdir}"
        ") 1>{log} 2>&1"


rule download_fastq_se:
    output:
        temp("data/se/{accession}.fastq.gz"),
    log:
        "results/logs/download_fastq/se/{accession}.log",
    params:
        out_prefix="data/se",
        extra=config["download_fastq_se"]["params"],
        tempdir=lambda wildcards: "data/se/{}.fasterq".format(wildcards.accession),
    threads: config["download_fastq_se"]["threads"]
    resources:
        max_downloads=1
    retries: config["download_fastq_se"]["retries"]
    conda:
        "../envs/sra-tools.yaml"
    shell:
        "("
        "fasterq-dump"
        " --skip-technical --force"
        " --threads {threads}"
        " --outdir {params.out_prefix}"
        " --temp {params.tempdir}"
        " {params.extra}"
        " {wildcards.accession}"
        "; pigz --processes {threads} {params.out_prefix}/{wildcards.accession}.fastq"
        "; rm -fr {params.tempdir}"
        ") 1>{log} 2>&1"


