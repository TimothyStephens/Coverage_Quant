

rule ref_parse_input:
    output:
        "resources/{ref_name}/ref.fa",
    log:
        "results/logs/resources/{ref_name}/ref_parse_genome.log",
    params:
        ref_file=lambda w: config["ref_path"][w.ref_name],
    conda:
        "../envs/bash.yaml"
    shell:
        "( if [[ {params.ref_file} == *.gz ]]; then zcat {params.ref_file} > {output}; else cat {params.ref_file} > {output}; fi ) 1>{log} 2>&1"


rule ref_faidx:
    input:
        "resources/{ref_name}/ref.fa",
    output:
        "resources/{ref_name}/ref.fa.fai",
    log:
        "results/logs/resources/{ref_name}/ref_faidx.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx"
        " {input}"
        " 1>{log} 2>&1"


rule ref_bbmap_index:
    input:
        "resources/{ref_name}/ref.fa",
    output:
        directory("resources/{ref_name}/ref.fa.idx")
    log:
        "results/logs/resources/{ref_name}/ref_bbmap_index.log",
    params:
        extra=config["ref_bbmap_index"]["params"],
    threads: config["ref_bbmap_index"]["threads"]
    resources:
        mem_gb=config["ref_bbmap_index"]["memory"]
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbmap.sh"
        " -Xmx{resources.mem_gb}g"
        " ref={input}"
        " path={output}"
        " overwrite=true"
        " build=1"
        " threads={threads}"
        " {params.extra}"
        " 1>{log} 2>&1"


