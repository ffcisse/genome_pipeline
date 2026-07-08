rule visuals:
    """
    Stage 6 -- Generate summary figures from the cross-genome summary tables
    (Phase 1: the boxplots/KDEs/heatmaps validated in the exploratory notebooks).

    Depends on `summaries` (already an aggregate over all genomes), so this rule
    stays a simple one-to-one dependency, no expand() needed here.
    """
    input:
        OUTDIR + "/summaries/summaries.done",
    output:
        OUTDIR + "/visuals/visuals.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
