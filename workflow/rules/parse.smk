rule parse:
    """
    Stage 2 -- Parse raw protein + CDS FASTA into canonical per-genome tables
    (real Bio.SeqIO parsing via workflow/scripts/parse.py, same approach as
    the exploratory notebooks: terminal stop codons stripped from protein
    sequences, CDS sequences uppercased, length recorded for both).

    Depends on `qc` for the same genome: Snakemake sees that this rule's `qc`
    input matches the `qc` rule's output pattern, and automatically schedules qc
    to run first for that genome.

    Outputs both the two canonical tables and a ".parsed.done" marker from
    the original Phase 0 scaffold. protein_properties/cds_properties read
    the real tables directly now, not the marker -- it's kept only because
    nothing has needed removing it, not because anything still depends on
    it.
    """
    input:
        qc=OUTDIR + "/qc/{genome}.qc.done",
        protein=lambda wildcards: PROTEIN_STAGED[wildcards.genome],
        cds=lambda wildcards: CDS_STAGED[wildcards.genome],
    output:
        protein_table=OUTDIR + "/parsed/{genome}/protein_table.csv",
        cds_table=OUTDIR + "/parsed/{genome}/cds_table.csv",
        done=OUTDIR + "/parsed/{genome}.parsed.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "python workflow/scripts/parse.py --genome {wildcards.genome} "
        "--protein {input.protein} --cds {input.cds} "
        "--protein-out {output.protein_table} --cds-out {output.cds_table} "
        "--done {output.done}"
