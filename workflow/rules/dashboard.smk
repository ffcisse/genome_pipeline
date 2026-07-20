rule dashboard_data:
    """
    Stage 7 (Phase 6a) -- Build the JSON data payload the dashboard HTML
    embeds inline. See workflow/scripts/build_dashboard_data.py's docstring
    for the exact-vs-sampled distinction this payload encodes.

    Reads every genome's qc.done report directly (expand(..., genome=GENOMES),
    same fan-in idiom as merge_protein_table's protein/disorder inputs) so
    the dashboard's Overview/QC numbers come from Phase 1's real per-genome
    output, not recomputed. Phase 4's effect_sizes_<grouping>.csv/
    sensitivity_leave_one_out.csv are passed through as-is (no stats
    recomputed here, same discipline as plot_effect_sizes.py/
    plot_sensitivity.py).

    params mirrors plot_boxplots.py's: primary/subgroup grouping column
    NAMES and their configured VALUE order come from config.yaml (via the
    Snakefile's SENSITIVITY_*/GROUP_VALUE_ORDER), nothing hardcoded here.
    DASHBOARD_SAMPLE_PER_GENOME/SEED are config.yaml's dashboard: block
    (Snakefile).
    """
    input:
        protein=OUTDIR + "/summaries/master_protein_table.csv",
        genomes_tsv=config["genome_table"],
        qc_reports=expand(OUTDIR + "/qc/{genome}.qc.done", genome=GENOMES),
        effect_sizes_primary=OUTDIR + "/summaries/effect_sizes_{}.csv".format(SENSITIVITY_PRIMARY_GROUPING),
        effect_sizes_subgroup=OUTDIR + "/summaries/effect_sizes_{}.csv".format(SENSITIVITY_SUBGROUP_COLUMN),
        sensitivity=OUTDIR + "/summaries/sensitivity_leave_one_out.csv",
    output:
        OUTDIR + "/dashboard/data.json",
    params:
        primary_grouping=SENSITIVITY_PRIMARY_GROUPING,
        subgroup_column=SENSITIVITY_SUBGROUP_COLUMN,
        primary_order=GROUP_VALUE_ORDER.get(SENSITIVITY_PRIMARY_GROUPING, []),
        subgroup_order=GROUP_VALUE_ORDER.get(SENSITIVITY_SUBGROUP_COLUMN, []),
        exclude_properties=VISUALS_EXCLUDE_PROPERTIES,
        sample_per_genome=DASHBOARD_SAMPLE_PER_GENOME,
        sample_seed=DASHBOARD_SAMPLE_SEED,
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/build_dashboard_data.py "
        "--master-protein-table {input.protein} --genomes-tsv {input.genomes_tsv} "
        "--qc-reports {input.qc_reports} "
        "--effect-sizes-primary {input.effect_sizes_primary} --effect-sizes-subgroup {input.effect_sizes_subgroup} "
        "--sensitivity-leave-one-out {input.sensitivity} "
        "--primary-grouping {params.primary_grouping} --subgroup-column {params.subgroup_column} "
        "--primary-order {params.primary_order} --subgroup-order {params.subgroup_order} "
        "--exclude-properties {params.exclude_properties} "
        "--sample-per-genome {params.sample_per_genome} --sample-seed {params.sample_seed} "
        "--output {output}"


rule dashboard:
    """
    Stage 7 (Phase 6a) -- Assemble the standalone dashboard HTML: inlines
    the vendored Plotly bundle + dashboard_data's JSON payload + this
    pipeline's custom HTML/CSS/JS (workflow/resources/dashboard/
    dashboard_template.html) into one file, so it opens directly in a
    browser with no server, no install, and no network access (see
    workflow/scripts/build_dashboard_html.py). Replaces the Phase 0 stub --
    this output path (results/dashboard/proteome_dashboard.html) is the
    real, final pipeline deliverable, so `rule all` in the Snakefile
    targets it directly rather than a separate touch-only marker.
    """
    input:
        data=OUTDIR + "/dashboard/data.json",
        plotly_js="workflow/resources/vendor/plotly-cartesian.min.js",
        template="workflow/resources/dashboard/dashboard_template.html",
        app_js="workflow/resources/dashboard/dashboard_app.js",
    output:
        OUTDIR + "/dashboard/proteome_dashboard.html",
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/build_dashboard_html.py "
        "--data-json {input.data} --plotly-js {input.plotly_js} --template {input.template} "
        "--app-js {input.app_js} --output {output}"
