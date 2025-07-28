process PREPARE_INPUTS {

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
    def meta_mapper_arg = meta_mapper ? "--meta-mapper $meta_mapper" : ''
    """
    bifrost-httr prepare-inputs \
        --meta-data $meta_data \
        --counts $counts \
        --config $bifrost_config \
        --output-dir bifrost_inputs \
        $meta_mapper_arg \
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
