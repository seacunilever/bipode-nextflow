process PREPARE_INPUTS {
    publishDir "${params.outdir}/", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

    input:
    path meta_data
    path meta_mapper
    path counts
    path substances_cell_types

    output:
    path "bifrost_inputs/*.json", emit: prepared_inputs

    script:
    def args = task.ext.args ?: ''
    """
    prepare_bifrost_inputs.py \
        --meta-data $meta_data \
        --meta-mapper $meta_mapper \
        --counts $counts \
        --substances-cell-types $substances_cell_types \
        --output-dir bifrost_inputs \
        $args
    """
}
