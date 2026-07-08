rule parse:
    """
    Stage 2 -- Parse raw protein + CDS FASTA into a canonical per-species
    intermediate table (Phase 1: real Bio.SeqIO parsing, same approach as the
    exploratory notebooks).

    Depends on `qc` for the same species: Snakemake sees that this rule's `qc`
    input matches the `qc` rule's output pattern, and automatically schedules qc
    to run first for that species.
    """
    input:
        qc=OUTDIR + "/qc/{species}.qc.done",
        protein=config["input"]["protein_dir"] + "/{species}.fasta",
        cds=config["input"]["cds_dir"] + "/{species}.fasta",
    output:
        OUTDIR + "/parsed/{species}.parsed.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
