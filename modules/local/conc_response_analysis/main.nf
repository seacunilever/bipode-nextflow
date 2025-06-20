process CONC_RESPONSE_ANALYSIS {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

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
    def probe_files_extract = probes.collect { "Data/" + it + ".pkl" }.join(" ")
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir Data Samples Fits

    tar -zxf $all_probe_file -C Data/ $probe_files

    # Compile Stan model if .stan file is provided
    if [[ "$model" == *.stan ]]; then
        echo "Compiling Stan model..."
        model_executable="${model.baseName}"
        bifrost-httr compile-model "$model"
    else
        model_executable="$model"
    fi

    bifrost-httr run-analysis \
        --data-files $probe_files_extract \
        --model-executable \$model_executable \
        --n-cores $task.cpus \
        $args

    sleep 5

    cd Fits && find . -name "*.pkl" -print0 | tar $args2 -cf - --null -T - | gzip $args3 > ../${prefix}.tar.gz
    """
}
