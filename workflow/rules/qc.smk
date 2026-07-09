rule qc:
    """
    Stage 1 -- Quality control on the raw input FASTA for one genome.

    Reports sequence counts, non-standard residues, internal stop codons,
    empty/duplicate records (protein) and non-multiple-of-3/non-ATGC/empty/
    duplicate records (CDS) -- same checks validated in the exploratory
    notebooks (workflow/scripts/qc.py). Hard-fails only on problems that make
    a file unusable (nothing parsed, duplicate IDs); a handful of X's or
    internal stops is normal in a real proteome and just gets reported.

    `{genome}` is a wildcard: Snakemake fills it in with each value from GENOMES
    (defined in the Snakefile) when it needs to produce a qc output for that genome.
    One rule definition, one job per genome -- that's the point of wildcards.

    Input is PROTEIN_STAGED/CDS_STAGED (built in the Snakefile), not a literal
    ".fasta" path, so this works whether or not stage_inputs is in the DAG.
    """
    input:
        protein=lambda wildcards: PROTEIN_STAGED[wildcards.genome],
        cds=lambda wildcards: CDS_STAGED[wildcards.genome],
    output:
        OUTDIR + "/qc/{genome}.qc.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "python workflow/scripts/qc.py --genome {wildcards.genome} "
        "--protein {input.protein} --cds {input.cds} --output {output}"
