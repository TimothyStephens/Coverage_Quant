

rule mapping_bbmap_pe:
    input:
        reads=rules.trimming_pe_merge.output.trimmed,
        ref=rules.ref_parse_input.output,
        idx=rules.ref_bbmap_index.output,
    output:
        "results/mapping/{ref_name}/pe/{sample}.coordSorted.bam",
    log:
        general="results/logs/mapping/{ref_name}/pe/{sample}.0_general.log",
        reformat="results/logs/mapping/{ref_name}/pe/{sample}.1_reformat.log",
        mapping="results/logs/mapping/{ref_name}/pe/{sample}.2_mapping.log",
        sorting="results/logs/mapping/{ref_name}/pe/{sample}.3_sorting.log",
    params:
        # Params
        reformat_params=config["mapping_bbmap_pe"]["reformat_params"],
        mapping_params=config["mapping_bbmap_pe"]["mapping_params"],
        sorting_params=config["mapping_bbmap_pe"]["sorting_params"],
        # Threads
        reformat_threads=config["mapping_bbmap_pe"]["reformat_threads"],
        mapping_threads=config["mapping_bbmap_pe"]["mapping_threads"],
        sorting_threads=config["mapping_bbmap_pe"]["sorting_threads"],
        # Memory
        reformat_memory=config["mapping_bbmap_pe"]["reformat_memory"],
        mapping_memory=config["mapping_bbmap_pe"]["mapping_memory"],
        sorting_memory=config["mapping_bbmap_pe"]["sorting_memory"],
    priority: 20
    threads: config["mapping_bbmap_pe"]["reformat_threads"] +
            config["mapping_bbmap_pe"]["mapping_threads"] +
            config["mapping_bbmap_pe"]["sorting_threads"]
    resources:
        mem_gb=config["mapping_bbmap_pe"]["reformat_memory"] +
                config["mapping_bbmap_pe"]["mapping_memory"] +
                (config["mapping_bbmap_pe"]["sorting_memory"] * config["mapping_bbmap_pe"]["sorting_threads"])
    conda:
        "../envs/bbmap.yaml"
    shell:
        "("
        "reformat.sh"
        "  {params.reformat_params}"
        "  in={input.reads[0]} in2={input.reads[1]}"
        "  out=stdout.fq"
        "  threads={params.reformat_threads}"
        "  -Xmx{params.reformat_memory}g"
        "  2>{log.reformat}"
        " | bbmap.sh"
        "  {params.mapping_params}"
        "  in=stdin.fq"
        "  interleaved=t"
        "  ref={input.ref}"
        "  path={input.idx}"
        "  outm=stdout.sam"
        "  threads={params.mapping_threads}"
        "  -Xmx{params.mapping_memory}g"
        "  2>{log.mapping}"
        " | samtools sort"
        " {params.sorting_params}"
        "  -@ {params.sorting_threads}"
        "  -m {params.sorting_memory}G"
        "  -T {output}.tmp"
        "  -o {output} -"
        "  2>{log.sorting}"
        ")"
        " 1>{log.general} 2>&1"


rule mapping_bbmap_se:
    input:
        reads=rules.trimming_se_merge.output.trimmed,
        ref=rules.ref_parse_input.output,
        idx=rules.ref_bbmap_index.output,
    output:
        "results/mapping/{ref_name}/se/{sample}.coordSorted.bam",
    log:
        general="results/logs/mapping/{ref_name}/se/{sample}.0_general.log",
        reformat="results/logs/mapping/{ref_name}/se/{sample}.1_reformat.log",
        mapping="results/logs/mapping/{ref_name}/se/{sample}.2_mapping.log",
        sorting="results/logs/mapping/{ref_name}/se/{sample}.3_sorting.log",
    params:
        # Params
        reformat_params=config["mapping_bbmap_se"]["reformat_params"],
        mapping_params=config["mapping_bbmap_se"]["mapping_params"],
        sorting_params=config["mapping_bbmap_se"]["sorting_params"],
        # Threads
        reformat_threads=config["mapping_bbmap_se"]["reformat_threads"],
        mapping_threads=config["mapping_bbmap_se"]["mapping_threads"],
        sorting_threads=config["mapping_bbmap_se"]["sorting_threads"],
        # Memory
        reformat_memory=config["mapping_bbmap_se"]["reformat_memory"],
        mapping_memory=config["mapping_bbmap_se"]["mapping_memory"],
        sorting_memory=config["mapping_bbmap_se"]["sorting_memory"],
    priority: 20
    threads: config["mapping_bbmap_se"]["reformat_threads"] +
            config["mapping_bbmap_se"]["mapping_threads"] +
            config["mapping_bbmap_se"]["sorting_threads"]
    resources:
        mem_gb=config["mapping_bbmap_se"]["reformat_memory"] +
                config["mapping_bbmap_se"]["mapping_memory"] +
                (config["mapping_bbmap_se"]["sorting_memory"] * config["mapping_bbmap_se"]["sorting_threads"])
    conda:
        "../envs/bbmap.yaml"
    shell:
        "("
        "reformat.sh"
        "  {params.reformat_params}"
        "  in={input.reads[0]}"
        "  out=stdout.fq"
        "  threads={params.reformat_threads}"
        "  -Xmx{params.reformat_memory}g"
        "  2>{log.reformat}"
        " | bbmap.sh"
        "  {params.mapping_params}"
        "  in=stdin.fq"
        "  ref={input.ref}"
        "  path={input.idx}"
        "  outm=stdout.sam"
        "  threads={params.mapping_threads}"
        "  -Xmx{params.mapping_memory}g"
        "  2>{log.mapping}"
        " | samtools sort"
        " {params.sorting_params}"
        "  -@ {params.sorting_threads}"
        "  -m {params.sorting_memory}G"
        "  -T {output}.tmp"
        "  -o {output} -"
        "  2>{log.sorting}"
        ")"
        " 1>{log.general} 2>&1"


rule mapping_coverm_coverage:
    input:
        ref=rules.ref_parse_input.output,
        g2s=config["g2s_path"],
        bams=get_bbmap_bams,
    output:
        "results/{project}/mapping/{ref_name}/coverage.tsv.gz",
    params:
        extra=config["mapping_coverm_coverage"]["params"],
        stats=config["mapping_coverm_coverage"]["stats"],
    log:
        "results/logs/{project}/mapping/{ref_name}/coverage.log",
    threads: config["mapping_coverm_coverage"]["threads"]
    conda:
        "../envs/coverm.yaml"
    shell:
        "("
        "coverm genome"
        "  {params.extra}"
        "  --output-format sparse --methods {params.stats}"
        "  --genome-definition {input.g2s}"
        "  -t {threads}"
        "  -b {input.bams}"
        " | gzip -c"
        " > {output}"
        ")"
        " 1>{log} 2>&1"


