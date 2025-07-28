process COMPILE_STAN_MODEL {
    tag "${model ?: 'default model'}"
    label 'process_single'

    stageInMode 'copy'

    input:
    path model

    output:
    path "compiled_model/*", emit: compiled_model
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def model_arg = model ? "--model $model" : ''
    """
    mkdir -p compiled_model && cd compiled_model
    bifrost-httr compile-model $model_arg $args
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p compiled_model
    touch compiled_model/stub_model

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """
}
