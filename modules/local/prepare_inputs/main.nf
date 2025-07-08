process PREPARE_INPUTS {

    conda "bifrost-httr=0.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/76/76e8817651482fe89237efe5d385050d40144519c9f0c9fc5b0f9ee506292428/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.3.1--b4c49de956618921' }"

    cpus { 2 * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    path meta_data
    path meta_mapper
    path counts
    path bifrost_config

    output:
    path "bifrost_inputs/*.json", emit: prepared_inputs
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bifrost-httr prepare-inputs \
        --meta-data $meta_data \
        --counts $counts \
        --config $bifrost_config \
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
