"""
Rules: miRNA-seq and UMI-based miRNA-seq (Parts 15, 15.2)
Container: mirna_env (09_mirna.yml / containers/mirna.sif)
Tools: miRDeep2 (Bowtie v1 alignment), umi_tools (Part 15.2)
Input: single-end FASTQ (miRNA reads are short, typically not paired-end)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

MIRNA_ADAPTER = config["mirna"]["adapter"]
MIRNA_MIN_LEN = config["mirna"]["min_length"]
MIRNA_MAX_LEN = config["mirna"]["max_length"]


rule umi_extract:
    """Extract UMIs before trimming (Part 15.2 only; skip if umi.enabled=false)."""
    input:
        r1 = get_r1,
    output:
        r1_umi = temp(f"{OUTDIR}/trimmed/{{sample}}_umi_extracted.fastq.gz"),
    log:   f"{LOGDIR}/umi/{{sample}}_extract.log"
    conda: "../../environments/09_mirna.yml"
    container: "file://containers/mirna.sif"
    resources: mem_mb=4000, runtime=30
    params:
        pattern = config["umi"]["pattern"],
    shell:
        """
        if [ "{config[umi][enabled]}" = "True" ]; then
            umi_tools extract \
                --stdin      {input.r1} \
                --stdout     {output.r1_umi} \
                --extract-method=string \
                --bc-pattern {params.pattern} \
                --log        {log}
        else
            cp {input.r1} {output.r1_umi}
        fi
        """


rule trim_mirna:
    """Trim adapter and filter by length for miRNA-seq."""
    input:
        r1 = f"{OUTDIR}/trimmed/{{sample}}_umi_extracted.fastq.gz",
    output:
        r1 = temp(f"{OUTDIR}/trimmed/{{sample}}_mirna_trimmed.fastq.gz"),
    log:   f"{LOGDIR}/trim_mirna/{{sample}}.log"
    conda: "../../environments/09_mirna.yml"
    container: "file://containers/mirna.sif"
    threads: 4
    resources: mem_mb=4000, runtime=30
    shell:
        """
        cutadapt \
            -a {MIRNA_ADAPTER} \
            --minimum-length {MIRNA_MIN_LEN} \
            --maximum-length {MIRNA_MAX_LEN} \
            --cores {threads} \
            -o {output.r1} \
            {input.r1} \
            2>{log}
        """


rule mirdeep2_mapper:
    """Collapse reads and align to genome with Bowtie v1 for miRDeep2."""
    input:
        r1    = f"{OUTDIR}/trimmed/{{sample}}_mirna_trimmed.fastq.gz",
        fasta = config["genome"]["fasta"],
    output:
        reads  = f"{OUTDIR}/mirna/{{sample}}_collapsed.fa",
        arf    = f"{OUTDIR}/mirna/{{sample}}_reads_vs_genome.arf",
        bowtie_index = temp(directory(f"{OUTDIR}/mirna/bowtie_idx_{{sample}}")),
    log:   f"{LOGDIR}/mirdeep2/{{sample}}_mapper.log"
    conda: "../../environments/09_mirna.yml"
    container: "file://containers/mirna.sif"
    threads: 8
    resources: mem_mb=16000, runtime=60
    shell:
        """
        # Build Bowtie v1 index from genome FASTA
        mkdir -p {output.bowtie_index}
        bowtie-build {input.fasta} {output.bowtie_index}/genome 2>>{log}

        # mapper.pl: gunzip, collapse, align (1 mismatch, -m 5 = discard >5 loci)
        mapper.pl \
            <(zcat {input.r1}) \
            -e -h -j -l {MIRNA_MIN_LEN} \
            -m -q \
            -p {output.bowtie_index}/genome \
            -s {output.reads} \
            -t {output.arf} \
            -o {threads} \
            2>>{log}
        """


rule mirdeep2_quantifier:
    """Run miRDeep2 core algorithm and quantify known miRNAs."""
    input:
        reads    = f"{OUTDIR}/mirna/{{sample}}_collapsed.fa",
        arf      = f"{OUTDIR}/mirna/{{sample}}_reads_vs_genome.arf",
        mature   = config["mirna"]["mirbase_mature"],
        hairpin  = config["mirna"]["mirbase_hairpin"],
        fasta    = config["genome"]["fasta"],
    output:
        html  = f"{OUTDIR}/mirna/{{sample}}_miRDeep2_results.html",
        tsv   = f"{OUTDIR}/mirna/{{sample}}_miRDeep2_results.tsv",
    log:   f"{LOGDIR}/mirdeep2/{{sample}}_quantifier.log"
    conda: "../../environments/09_mirna.yml"
    container: "file://containers/mirna.sif"
    threads: 4
    resources: mem_mb=16000, runtime=90
    params:
        outdir  = f"{OUTDIR}/mirna",
        species = lambda wc: "hsa" if config["organism"] == "human" else "mmu",
    shell:
        """
        cd {params.outdir}
        miRDeep2.pl \
            {input.reads} \
            {input.fasta} \
            {input.arf} \
            {input.mature} \
            none \
            {input.hairpin} \
            -P -v \
            2>{log}
        # miRDeep2 names output by timestamp; rename to deterministic names
        mv result_*.html {wildcards.sample}_miRDeep2_results.html 2>>{log} || true
        mv result_*.tsv  {wildcards.sample}_miRDeep2_results.tsv  2>>{log} || true
        """
