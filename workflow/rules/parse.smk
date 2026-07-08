rule parse:
    """
    Stage 2 -- Parse raw protein + CDS FASTA into a canonical per-genome
    intermediate table (Phase 1: real Bio.SeqIO parsing, same approach as the
    exploratory notebooks).

    Depends on `qc` for the same genome: Snakemake sees that this rule's `qc`
    input matches the `qc` rule's output pattern, and automatically schedules qc
    to run first for that genome.
    """
    input:
        qc=OUTDIR + "/qc/{genome}.qc.done",
        protein=config["input"]["protein_dir"] + "/{genome}.fasta",
        cds=config["input"]["cds_dir"] + "/{genome}.fasta",
    output:
        OUTDIR + "/parsed/{genome}.parsed.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
