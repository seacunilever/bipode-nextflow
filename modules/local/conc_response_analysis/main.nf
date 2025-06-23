process CONC_RESPONSE_ANALYSIS {
    tag "${meta.id}"

    conda "bifrost-httr=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05fa08012fb11ccb282c05e1b48b53c6220b0853692cafcbaf8829749d6aabc/data' :
        'wave.seqera.io/wt/25fa77f460cd/wave/build:bifrost-httr-0.1.0--2c648d2de87966a9' }"

    input:
    path model
    tuple val(meta), val(probes), path(all_probe_file)

    output:
    tuple val(meta), path("${prefix}.tar.gz"), emit: compressed_fits_files
    tuple val(meta), path("Fits/*.pkl"), emit: fits_files
    tuple val(meta), path("Fits/*.json"), emit: json_summaries

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def probe_files = probes.collect { it + ".pkl" }.join(" ")
    def probe_files_extract = probes.collect { " -f Data/" + it + ".pkl" }.join(" ")
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir Data Samples Fits

    tar -zxf $all_probe_file -C Data/ $probe_files

    # Handle model: empty (compile default), .stan (compile file), or executable (use as-is)
    if [[ -z "$model" ]]; then
        mkdir -p ModelCompilation && cd ModelCompilation
        bifrost-httr compile-model
        executable_file=\$(find . -maxdepth 1 -type f -perm +111 | head -n 1)
        if [[ -z "\$executable_file" ]]; then
            echo "Error: No executable found after compilation"
            exit 1
        fi
        model_executable=\$(realpath "\$executable_file")
        cd ..
    elif [[ "$model" == *.stan ]]; then
        model_executable="${model.baseName}"
        bifrost-httr compile-model "$model"
    else
        model_executable="$model"
    fi

    bifrost-httr run-analysis \
        $probe_files_extract \
        --model-executable \$model_executable \
        --n-cores $task.cpus \
        $args

    sleep 5

    cd Fits && find . -name "*.pkl" -print0 | tar $args2 -cf - --null -T - | gzip $args3 > ../${prefix}.tar.gz
    """
}
