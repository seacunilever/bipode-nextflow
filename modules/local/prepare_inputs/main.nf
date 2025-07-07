process PREPARE_INPUTS {

    conda "bifrost-httr=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37b8df53325072038a7dce735a434031e56946aeca172a683216dbdbc197fc5d/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.3.0--6161c4cc71c68c4c' }"

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
        --counts $counts \
        --config $substances_cell_types \
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
