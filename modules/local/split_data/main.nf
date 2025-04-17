process SPLIT_DATA {
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

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
