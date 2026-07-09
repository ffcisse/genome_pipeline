rule stage_inputs:
    """
    Stage 0 -- Symlink each genome's real protein + CDS FASTA (wherever it
    actually lives -- e.g. several directories deep in a Mycocosm download
    tree, gzipped or gzipped-tarred) into a flat layout under
    input.protein_dir / input.cds_dir, named "<genome_id><real extension>".

    Only included when config["staging"]["source_dir"] is set (see
    config.yaml). Symlinks, not copies: no data duplication, and re-running
    after the source data changes on disk just re-points the link.

    PROTEIN_SOURCE/CDS_SOURCE (real files) and PROTEIN_STAGED/CDS_STAGED
    (symlink destinations, with the real extension already resolved) are
    built once in the Snakefile, before this rule is included, since the
    source files already exist on disk -- no checkpoint needed to learn their
    extensions.
    """
    input:
        protein=list(PROTEIN_SOURCE.values()),
        cds=list(CDS_SOURCE.values()),
    output:
        protein=list(PROTEIN_STAGED.values()),
        cds=list(CDS_STAGED.values()),
    run:
        for genome in GENOMES:
            for src, dst in (
                (PROTEIN_SOURCE[genome], PROTEIN_STAGED[genome]),
                (CDS_SOURCE[genome], CDS_STAGED[genome]),
            ):
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                src_abs = os.path.abspath(src)
                if os.path.islink(dst):
                    if os.path.realpath(dst) == os.path.realpath(src_abs):
                        continue
                    os.remove(dst)
                elif os.path.exists(dst):
                    raise FileExistsError(
                        f"{dst} already exists and is not a symlink -- refusing to overwrite"
                    )
                os.symlink(src_abs, dst)
