rule dashboard:
    """
    Stage 7 -- Build the interactive dashboard.

    Deferred: this is intentionally just an empty rule + output path for now
    (per scaffold scope). No conda env declared yet since there's no real
    dependency to pin until dashboard internals are designed.
    """
    input:
        OUTDIR + "/visuals/visuals.done",
    output:
        OUTDIR + "/dashboard/dashboard.done",
    shell:
        "touch {output}"
