process SPLIT_DATA {
    input:
    tuple val(name), path(input_data)

    output:
    tuple val(name), path("${name}_probes.tar.gz"), emit: all_probe_files

    script:
    """
    split_data.py --input-file $input_data --analysis-dir .
    tar -czf ${name}_probes.tar.gz -C Data/ .
    """
} 