rule merge_protein_table:
    """
    Stage 5 (Phase 4) -- Join every genome's protein_properties.csv +
    disorder.csv into one master table, with group labels joined in from
    config/genomes.tsv (which is why it's listed as an input here too, not
    just read incidentally -- editing it should invalidate this rule). See
    workflow/scripts/merge_protein_table.py's docstring for the row-count
    validation (must exactly match the concatenated input, checked via both
    pandas' `validate="one_to_one"` and an explicit count comparison --
    fails loudly, not a warning, on any mismatch).

    params.group_columns (GROUPING_COLUMNS, from the Snakefile) is exactly
    config.yaml's sensitivity.primary_grouping/subgroup_column -- which
    genome_table columns actually get attached to the master table. Nothing
    here hardcodes "lifestyle"/"lineage".

    This is a FAN-IN rule like the old `summaries` stub was: expand(...)
    turns the {genome} wildcard template into one concrete path per genome
    in GENOMES, so Snakemake won't run this until every genome has finished
    both protein_properties and disorder.
    """
    input:
        protein=expand(OUTDIR + "/protein_properties/{genome}/protein_properties.csv", genome=GENOMES),
        disorder=expand(OUTDIR + "/disorder/{genome}/disorder.csv", genome=GENOMES),
        genomes_tsv=config["genome_table"],
    output:
        OUTDIR + "/summaries/master_protein_table.csv",
    params:
        group_columns=GROUPING_COLUMNS,
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/merge_protein_table.py "
        "--protein-tables {input.protein} --disorder-tables {input.disorder} "
        "--genomes-tsv {input.genomes_tsv} --group-columns {params.group_columns} --output {output}"


rule merge_cds_table:
    """
    Stage 5 (Phase 4) -- Same as merge_protein_table, for the CDS side:
    joins cds_properties.csv + codon_usage.csv per genome, plus group
    labels (params.group_columns, same as merge_protein_table above). See
    workflow/scripts/merge_cds_table.py.
    """
    input:
        cds=expand(OUTDIR + "/cds_properties/{genome}/cds_properties.csv", genome=GENOMES),
        codon_usage=expand(OUTDIR + "/cds_properties/{genome}/codon_usage.csv", genome=GENOMES),
        genomes_tsv=config["genome_table"],
    output:
        OUTDIR + "/summaries/master_cds_table.csv",
    params:
        group_columns=GROUPING_COLUMNS,
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/merge_cds_table.py "
        "--cds-tables {input.cds} --codon-usage-tables {input.codon_usage} "
        "--genomes-tsv {input.genomes_tsv} --group-columns {params.group_columns} --output {output}"


rule species_summary:
    """
    Stage 5 (Phase 4) -- One row per genome: median/mean/std of every
    numeric property from both master tables, plus gene counts and group
    labels. See workflow/scripts/species_summary.py.
    """
    input:
        protein=OUTDIR + "/summaries/master_protein_table.csv",
        cds=OUTDIR + "/summaries/master_cds_table.csv",
        genomes_tsv=config["genome_table"],
    output:
        OUTDIR + "/summaries/species_summary.csv",
    params:
        group_columns=GROUPING_COLUMNS,
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/species_summary.py "
        "--master-protein-table {input.protein} --master-cds-table {input.cds} "
        "--genomes-tsv {input.genomes_tsv} --group-columns {params.group_columns} --output {output}"


rule effect_sizes:
    """
    Stage 5 (Phase 4) -- Group-comparison effect sizes (Mann-Whitney U,
    CLES, rank-biserial correlation) for every numeric property, for a
    CONFIGURABLE grouping variable -- this is the key requirement Phase 4
    was built around. See workflow/scripts/effect_sizes.py's docstring for
    why: with n~61,000, p-values are ~0 regardless of whether a difference
    is biologically meaningful, and some properties turn out to separate by
    a finer subgroup rather than the primary grouping -- so being able to
    swap the grouping variable without touching code is the whole point.
    Nothing here or in effect_sizes.py hardcodes which columns or values
    those are (an earlier version did -- see git history); group NAMES come
    from GROUPING_COLUMNS (config.yaml's sensitivity: block), and group
    VALUE order comes from params.group_order (config.yaml's optional
    group_value_order: block, empty/alphabetical if that column isn't
    listed there).

    `{grouping}` is a Snakemake wildcard, same mechanism as `{genome}`
    elsewhere in this pipeline, but instead of matching a value from
    GENOMES, it's constrained below to GROUPING_COLUMNS via
    `wildcard_constraints` (a regex Snakemake requires the wildcard's value
    to match) -- built from that config-derived list, not a hardcoded
    "lifestyle|lineage" string, so it automatically reflects whatever
    config.yaml's sensitivity block actually names. Without this constraint,
    a typo'd target like `results/summaries/effect_sizes_lifestyle_typo.csv`
    would still match this rule at DAG-build time (Snakemake would only
    catch the problem later, at runtime, when effect_sizes.py tried to find
    a "lifestyle_typo" column) -- the constraint catches it immediately
    instead, with a clearer error.

    Requesting both variants: expand(... "{grouping}" ..., grouping=[...])
    in the `summaries` rule below is exactly the same expand() mechanism as
    GENOMES, just over GROUPING_COLUMNS (a 2-item list by default) instead
    of the 9 genomes -- it's what makes both effect_sizes_lifestyle.csv and
    effect_sizes_lineage.csv part of the default build.
    """
    input:
        protein=OUTDIR + "/summaries/master_protein_table.csv",
        cds=OUTDIR + "/summaries/master_cds_table.csv",
        genomes_tsv=config["genome_table"],
    output:
        OUTDIR + "/summaries/effect_sizes_{grouping}.csv",
    wildcard_constraints:
        grouping="|".join(re.escape(g) for g in GROUPING_COLUMNS),
    params:
        group_order=lambda wildcards: GROUP_VALUE_ORDER.get(wildcards.grouping, []),
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/effect_sizes.py --grouping {wildcards.grouping} "
        "--genomes-tsv {input.genomes_tsv} --group-order {params.group_order} "
        "--master-protein-table {input.protein} --master-cds-table {input.cds} --output {output}"


rule sensitivity_drop_subgroup:
    """
    Stage 5 (Phase 4) -- Leave-one-subgroup-out sensitivity analysis:
    quantifies how much of each property's SENSITIVITY_PRIMARY_GROUPING
    effect size is actually driven by one SENSITIVITY_SUBGROUP_COLUMN value
    (both configured in config.yaml's `sensitivity:` block, e.g. by default
    "is the apparent lifestyle effect actually a specific lineage's
    artifact?"). See workflow/scripts/sensitivity_drop_subgroup.py --
    nothing here or in that script is hardcoded to any one dataset's group
    names or values.

    `{subgroup}` is a wildcard, same mechanism as `{genome}`/`{grouping}`
    elsewhere in this pipeline. Its wildcard_constraints is built FROM THE
    DATA (SUBGROUPS, computed in the Snakefile from genome_table's
    SENSITIVITY_SUBGROUP_COLUMN column) rather than a hardcoded regex like
    effect_sizes' -- so it automatically reflects whatever subgroup values
    genome_table actually has for this dataset, with no edit needed here if
    genomes are added/removed/relabeled.

    `params:` (as opposed to `input:`) is how a rule passes plain
    configuration values -- not files Snakemake needs to track for the DAG
    -- into its shell command. primary_grouping/subgroup_column/group_order
    are the same for every {subgroup} value, so they're fixed params here
    rather than wildcards.
    """
    input:
        protein=OUTDIR + "/summaries/master_protein_table.csv",
        cds=OUTDIR + "/summaries/master_cds_table.csv",
        genomes_tsv=config["genome_table"],
    output:
        OUTDIR + "/summaries/sensitivity_drop_{subgroup}.csv",
    wildcard_constraints:
        subgroup="|".join(re.escape(s) for s in SUBGROUPS),
    params:
        primary_grouping=SENSITIVITY_PRIMARY_GROUPING,
        subgroup_column=SENSITIVITY_SUBGROUP_COLUMN,
        group_order=GROUP_VALUE_ORDER.get(SENSITIVITY_PRIMARY_GROUPING, []),
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/sensitivity_drop_subgroup.py "
        "--primary-grouping {params.primary_grouping} --subgroup-column {params.subgroup_column} "
        "--exclude-subgroup {wildcards.subgroup} --genomes-tsv {input.genomes_tsv} "
        "--group-order {params.group_order} "
        "--master-protein-table {input.protein} --master-cds-table {input.cds} --output {output}"


rule sensitivity_leave_one_out:
    """
    Stage 5 (Phase 4) -- Combines every subgroup's sensitivity_drop_{subgroup}.csv
    into one sweep table (one row per property per excluded subgroup), so
    you can see how much EACH subgroup individually contributes to the
    apparent primary-grouping effect, not just one hand-picked subgroup.

    input: uses expand(..., subgroup=SUBGROUPS) to depend on every
    individual sensitivity_drop_subgroup output -- the same fan-in idiom as
    merge_protein_table's expand(..., genome=GENOMES), just over SUBGROUPS
    instead. params.subgroups passes that same list (in the same order) to
    the combining script, so it can label each file's rows with the
    subgroup that produced them without having to parse it back out of the
    filename.
    """
    input:
        expand(OUTDIR + "/summaries/sensitivity_drop_{subgroup}.csv", subgroup=SUBGROUPS),
    output:
        OUTDIR + "/summaries/sensitivity_leave_one_out.csv",
    params:
        subgroups=SUBGROUPS,
    conda:
        "../envs/summaries.yaml"
    shell:
        "python workflow/scripts/combine_sensitivity_sweep.py "
        "--inputs {input} --subgroups {params.subgroups} --output {output}"


rule summaries:
    """
    Stage 5 checkpoint -- Phase 4 is summaries ONLY (visuals/dashboard come
    later), so this stays a thin touch-only aggregator, exactly like it was
    as a Phase 0 stub: visuals.smk/dashboard.smk are untouched by Phase 4,
    they still just need results/summaries/summaries.done to exist. Its
    `input:` is what actually pulls in every real Phase 4 rule above (and,
    transitively through THOSE rules' own inputs, protein_properties/
    disorder/cds_properties/codon_usage for every genome) -- the touch
    itself doesn't do anything, but Snakemake won't run it (or let anything
    downstream proceed) until every input listed here exists.

    effect_sizes' expand() below uses GROUPING_COLUMNS (config-derived), not
    a hardcoded ["lifestyle", "lineage"] list -- so the default build
    always requests exactly the two effect_sizes tables that make sense for
    whatever config.yaml's sensitivity block names, whatever this
    deployment's genome_table calls them.
    """
    input:
        master_protein=OUTDIR + "/summaries/master_protein_table.csv",
        master_cds=OUTDIR + "/summaries/master_cds_table.csv",
        species_summary=OUTDIR + "/summaries/species_summary.csv",
        effect_sizes=expand(OUTDIR + "/summaries/effect_sizes_{grouping}.csv", grouping=GROUPING_COLUMNS),
        sensitivity=OUTDIR + "/summaries/sensitivity_leave_one_out.csv",
    output:
        OUTDIR + "/summaries/summaries.done",
    shell:
        "touch {output}"
