process CONC_RESPONSE_ANALYSIS {
    tag "${meta.id}"

    conda "bifrost-httr=0.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed90af4777d8d7086ed99d0a825f99e20e39278a75c70d3f4f7b6336edf7e210/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.2.0--e8ca5c015e9a6142' }"

    cpus { params.n_cores }
    memory { 3.GB * task.attempt }
    time { 4.h * task.attempt }
    disk '10 GB'

    input:
    path model
    tuple val(meta), val(probes), path(all_probe_file)

    output:
    tuple val(meta), path("${prefix}.tar.gz"), emit: compressed_fits_files
    tuple val(meta), path("Fits/*.pkl"), emit: fits_files
    tuple val(meta), path("Fits/*.json"), emit: json_summaries
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def probe_files = probes.collect { it + ".pkl" }.join(" ")
    def probe_files_extract = probes.collect { " -f Data/" + it + ".pkl" }.join(" ")
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p Data Samples Fits

    tar -zxf $all_probe_file -C Data/ $probe_files

    # Handle model: empty (compile default), .stan (compile file), or executable (use as-is)
    if [[ -z "$model" ]] || [[ "$model" == *.stan ]]; then
        model_input=\${model:-""}  # Use empty string if model is empty, otherwise use model path
        bifrost-httr compile-model \$model_input
        model_executable=\$(find . -maxdepth 1 -type f -exec test -x {} \\; -print | head -n 1)
        if [[ -z "\$model_executable" ]]; then
            echo "Error: No executable found after compilation"
            exit 1
        fi
    else
        model_executable="$model"
    fi

    bifrost-httr run-analysis \
        $probe_files_extract \
        --model-executable \$model_executable \
        --n-cores $task.cpus \
        $args

    sleep 5

    cd Fits && find . -name "*.pkl" -print0 | tar $args2 -cf - --null -T - | gzip $args3 > ../${prefix}.tar.gz && cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p Fits
    touch Fits/stub.pkl
    touch Fits/stub.json
    touch ${prefix}.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """
}
