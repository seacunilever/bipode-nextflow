process PREPARE_INPUTS {

    conda "bifrost-httr=0.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed90af4777d8d7086ed99d0a825f99e20e39278a75c70d3f4f7b6336edf7e210/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.2.0--e8ca5c015e9a6142' }"

    cpus { 2 * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    path meta_data
    path meta_mapper
    path counts
    path substances_cell_types

    output:
    path "bifrost_inputs/*.json", emit: prepared_inputs
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bifrost-httr prepare-inputs \
        --meta-data $meta_data \
        --meta-mapper $meta_mapper \
        --counts $counts \
        --substances-cell-types $substances_cell_types \
        --output-dir bifrost_inputs \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bifrost_inputs
    touch bifrost_inputs/stub.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """
}
