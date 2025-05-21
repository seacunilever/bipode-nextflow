process CREATE_REPORTS {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "bifrost-reporting"

    input:
    tuple val(meta), path(input_file), val(cell_type), val(test_substance)

    output:
    path "${prefix}.pdf", emit: report_pdf

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_reports.py \
        --summary-file \"$input_file\" \
        --output-name \"$prefix\" \
        --test-substance \"$test_substance\" \
        --cell-type \"$cell_type\" \
        $args
    """
}
