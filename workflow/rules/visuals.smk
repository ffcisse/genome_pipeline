rule plot_boxplots:
    """
    Stage 6 (Phase 5) -- Boxplots + violin plots of every protein property:
    by primary_grouping, by subgroup_column, and by species (colored by
    subgroup, median-labeled). See workflow/scripts/plot_boxplots.py --
    genome-agnostic by construction, same discipline as Phase 4: property
    list, group values/order/colors, and species order all come from the
    data + config.yaml (params below), nothing hardcoded to this
    deployment's actual property/genome/group names.

    output: directory(...) rather than individual filenames -- the actual
    number of figures depends on how many numeric properties the master
    table happens to have, which isn't knowable at DAG-build time without
    a Snakemake `checkpoint` (a heavier mechanism for when the DAG itself
    must branch on a rule's output content; not needed here, since every
    plot rule's *input* is already fully known -- only the *count of
    output files* is data-dependent). Snakemake tracks a directory() output
    by its own existence/mtime, the same idea as this pipeline's .done
    marker files, just for a rule whose real product is "a batch of files
    in one place" rather than one countable file.

    params (not input:) is how config-driven values that aren't files --
    grouping column names, their configured value order, which properties
    to skip -- reach the script, same as Phase 4's rules.
    """
    input:
        protein=OUTDIR + "/summaries/master_protein_table.csv",
        genomes_tsv=config["genome_table"],
    output:
        directory(OUTDIR + "/plots/boxplots"),
    params:
        primary_grouping=SENSITIVITY_PRIMARY_GROUPING,
        subgroup_column=SENSITIVITY_SUBGROUP_COLUMN,
        primary_order=GROUP_VALUE_ORDER.get(SENSITIVITY_PRIMARY_GROUPING, []),
        subgroup_order=GROUP_VALUE_ORDER.get(SENSITIVITY_SUBGROUP_COLUMN, []),
        exclude_properties=VISUALS_EXCLUDE_PROPERTIES,
    conda:
        "../envs/visuals.yaml"
    shell:
        "python workflow/scripts/plot_boxplots.py "
        "--master-protein-table {input.protein} --genomes-tsv {input.genomes_tsv} "
        "--primary-grouping {params.primary_grouping} --subgroup-column {params.subgroup_column} "
        "--primary-order {params.primary_order} --subgroup-order {params.subgroup_order} "
        "--exclude-properties {params.exclude_properties} "
        "--output-dir {output}"


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
