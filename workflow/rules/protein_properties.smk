rule protein_properties:
    """
    Stage 3 (Phase 2a) -- Compute lightweight per-protein physicochemical
    properties for one genome: amino acid composition, pI, GRAVY, aliphatic
    index, thermostable fraction, instability index, net charge, charge
    density, cysteine fraction, carbon oxidation state, and AGGRESCAN-style
    aggregation propensity. All pure sequence/composition math -- fast,
    single-core, safe to run interactively. Ported from
    resources/01_proteome_overview_FINAL.ipynb -- see
    workflow/scripts/protein_properties.py.

    Intrinsic disorder prediction is deliberately NOT here -- it's heavy
    (needs its own SLURM resources) and lands in Phase 2b as a separate rule
    that reads this rule's output and adds a disorder column alongside it.

    Input is the real parsed table (not just parse's ".parsed.done" marker),
    since this rule actually needs the sequences, not just a signal that
    parsing finished. Snakemake still infers "run parse for this genome
    first" from the fact that this path is parse's declared output.
    """
    input:
        OUTDIR + "/parsed/{genome}/protein_table.csv",
    output:
        OUTDIR + "/protein_properties/{genome}/protein_properties.csv",
    conda:
        "../envs/analysis.yaml"
    shell:
        "python workflow/scripts/protein_properties.py --genome {wildcards.genome} "
        "--protein-table {input} --output {output}"
