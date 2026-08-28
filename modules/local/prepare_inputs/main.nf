process PREPARE_INPUTS {

    cpus 2
    memory 1.GB

    input:
    path meta_data
    path meta_mapper
    path counts
    path bipode_config

    output:
    path "bipode_inputs/*.json", emit: prepared_inputs
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def meta_mapper_arg = meta_mapper ? "--meta-mapper $meta_mapper" : ''
    """
    bipode-httr prepare-inputs \
        --meta-data $meta_data \
        --counts $counts \
        --config $bipode_config \
        --output-dir bipode_inputs \
        $meta_mapper_arg \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bipode-httr: \$(bipode-httr --version | sed 's/bipode-httr, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bipode_inputs
    touch bipode_inputs/stub.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bipode-httr: \$(bipode-httr --version | sed 's/bipode-httr, version //')
    END_VERSIONS
    """
}
