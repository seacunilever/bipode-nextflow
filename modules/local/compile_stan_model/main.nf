process COMPILE_STAN_MODEL {
    tag "$model"
    label 'process_single'

    conda "bifrost-httr=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05fa08012fb11ccb282c05e1b48b53c6220b0853692cafcbaf8829749d6aabc/data' :
        'wave.seqera.io/wt/25fa77f460cd/wave/build:bifrost-httr-0.1.0--2c648d2de87966a9' }"

    stageInMode 'copy'

    input:
    path model

    output:
    path "compiled_model/*", emit: compiled_model

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
    """
}
