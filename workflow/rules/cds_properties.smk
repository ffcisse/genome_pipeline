rule cds_properties:
    """
    Stage 4 (Phase 3) -- Compute per-gene CDS/codon-usage properties for one
    genome via codonW (Peden), the standard tool for ENC/GC/GC3s and codon
    usage in published papers -- NOT reimplemented in Python here. See
    workflow/scripts/cds_properties.py's docstring for exactly how its
    output is invoked and parsed (including why a synthetic FASTA header is
    used and why that's still a safe join back to the real cds_id).

    Depends on `parse` for the same genome. Runs independently of (and in
    parallel with) `protein_properties`/`disorder` -- all three only need
    `parse`'s output. Same idiom as protein_properties.smk: input is the
    real parsed table, not just parse's ".parsed.done" marker, since this
    rule actually needs the sequences; Snakemake still infers "run parse for
    this genome first" from the fact that this path is parse's declared
    output.

    Two outputs, not one, for the same reason Phase 2b's disorder rule
    writes a separate disorder.csv instead of extending protein_properties.csv:
    cds_properties.csv (one row per CDS: length, ENC, GC, GC3s, start/stop
    codon) stays narrow and readable, while the 64-column per-codon usage
    table lives in its own codon_usage.csv, joinable back on
    (genome, cds_id). Light/fast, like protein_properties -- codon math on
    parsed sequences, no heavy model like Phase 2b's disorder rule -- so no
    special `resources:` block is needed here.
    """
    input:
        OUTDIR + "/parsed/{genome}/cds_table.csv",
    output:
        properties=OUTDIR + "/cds_properties/{genome}/cds_properties.csv",
        codon_usage=OUTDIR + "/cds_properties/{genome}/codon_usage.csv",
    conda:
        "../envs/codonw.yaml"
    shell:
        "python workflow/scripts/cds_properties.py --genome {wildcards.genome} "
        "--cds-table {input} --properties-output {output.properties} "
        "--codon-usage-output {output.codon_usage}"
