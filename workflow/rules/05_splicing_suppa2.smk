"""
Rules: Differential alternative splicing with SUPPA2 (Part 9)
Environment: suppa2_env (06_suppa2.yml)
Tools: Salmon (transcript quant) → SUPPA2 PSI calculation → dpsi/dsigma test
Input: Salmon quant.sf files from rule salmon_quant (01_qc_align.smk)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


rule suppa2_generate_events:
    input:
        gtf = config["genome"]["gtf"],
    output:
        ioe = directory(f"{OUTDIR}/splicing/suppa2/events"),
    log:   f"{LOGDIR}/suppa2/generate_events.log"
    resources: mem_mb=8000, runtime=30
    params:
        activate = mamba_activate("suppa2_env"),
        outdir   = f"{OUTDIR}/splicing/suppa2/events",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {output.ioe}
        python -m suppa generateEvents \
            -i {input.gtf} \
            -o {params.outdir}/all_events \
            -f ioe \
            -e SE SS MX RI FL \
            2>{log}
        """


rule suppa2_psi_per_sample:
    input:
        quant = f"{OUTDIR}/salmon/{{sample}}/quant.sf",
        ioe   = f"{OUTDIR}/splicing/suppa2/events",
    output:
        psi = f"{OUTDIR}/splicing/suppa2/psi/{{sample}}.psi",
    log:   f"{LOGDIR}/suppa2/psi_{{sample}}.log"
    resources: mem_mb=8000, runtime=30
    params:
        activate = mamba_activate("suppa2_env"),
        events   = f"{OUTDIR}/splicing/suppa2/events/all_events",
    shell:
        """
        set -eo pipefail
        {params.activate}
        awk 'NR>1{{print $1"\t"$4}}' {input.quant} > /tmp/{wildcards.sample}_tpm.tsv
        python -m suppa psiPerEvent \
            -i {params.events}_*.ioe \
            -e /tmp/{wildcards.sample}_tpm.tsv \
            -o {output.psi} \
            2>{log}
        """


rule suppa2_dpsi:
    input:
        psi_all = expand(f"{OUTDIR}/splicing/suppa2/psi/{{sample}}.psi", sample=SAMPLES),
        ioe     = f"{OUTDIR}/splicing/suppa2/events",
    output:
        dpsi = f"{OUTDIR}/splicing/suppa2/{{comp}}_dpsi.tab",
        sig  = f"{OUTDIR}/splicing/suppa2/{{comp}}_significance.tab",
    log:   f"{LOGDIR}/suppa2/dpsi_{{comp}}.log"
    threads: 4
    resources: mem_mb=16000, runtime=60
    params:
        activate  = mamba_activate("suppa2_env"),
        events    = f"{OUTDIR}/splicing/suppa2/events/all_events",
        treatment = lambda wc: wc.comp.split("_vs_")[0],
        control   = lambda wc: wc.comp.split("_vs_")[1],
        outpfx    = f"{OUTDIR}/splicing/suppa2/{{comp}}",
    shell:
        """
        set -eo pipefail
        {params.activate}
        python -m suppa diffSplice \
            -m empirical \
            -i {params.events}_*.ioe \
            -p {params.outpfx} \
            -s {params.treatment} {params.control} \
            -gc \
            -o {params.outpfx} \
            2>{log}
        """
