process COMPILE_STAN_MODEL {
    tag "${model ?: 'default model'}"
    label 'process_single'

    conda "bifrost-httr=0.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed90af4777d8d7086ed99d0a825f99e20e39278a75c70d3f4f7b6336edf7e210/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.2.0--e8ca5c015e9a6142' }"

    stageInMode 'copy'

    input:
    path model

    output:
    path "compiled_model/*", emit: compiled_model
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p compiled_model && cd compiled_model

    # Handle model: empty (compile default) or .stan (compile file)
    if [[ -z "$model" ]]; then
        bifrost-httr compile-model $args
    else
        bifrost-httr compile-model "$model" $args
    fi
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
