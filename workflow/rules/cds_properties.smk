rule cds_properties:
    """
    Stage 4 -- Compute per-gene CDS properties for one species (GC, GC3, codon
    usage, ENC -- ported from the exploratory notebook in Phase 1, with codonW
    established there as the authoritative source for ENC/GC3/GC).

    Depends on `parse` for the same species. Runs independently of (and in
    parallel with) `protein_properties` -- both only need `parse`'s output.
    """
    input:
        OUTDIR + "/parsed/{species}.parsed.done",
    output:
        OUTDIR + "/cds_properties/{species}.properties.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
