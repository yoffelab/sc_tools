# sc_tools repo-root Snakefile — Lint, format, test
# Run from repo root: snakemake -s Snakefile [lint|format|test|all]

rule lint:
    output:
        touch(".snakemake_lint_done")
    shell:
        "ruff check sc_tools && ruff format --check sc_tools"

rule format:
    output:
        touch(".snakemake_format_done")
    shell:
        "ruff format sc_tools"

rule test:
    output:
        touch(".snakemake_test_done")
    shell:
        "python -m pytest sc_tools/tests -v"

rule all:
    input:
        ".snakemake_lint_done",
        ".snakemake_test_done"
    run:
        pass
