process CONC_RESPONSE_ANALYSIS {
    tag "${meta.id}"

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
    mkdir -p Data

    tar -zxf $all_probe_file -C Data/ $probe_files

    # Handle model: empty (compile default), .stan (compile file), or executable (use as-is)
    if [[ -z "$model" ]] || [[ "$model" == *.stan ]]; then
        model_input=\${model:-""}  # Use empty string if model is empty, otherwise use model path
        compilation_output=\$(bifrost-httr compile-model \$model_input 2>&1 | tee /dev/stderr)
        model_executable=\$(echo "\$compilation_output" | grep "Model compiled successfully:" | sed 's/.*Model compiled successfully: //')
        if [[ -z "\$model_executable" ]]; then
            echo "Error: No executable path found in compilation output"
            exit 1
        fi
    else
        model_executable="$model"
    fi

    bifrost-httr run-analysis \
        $probe_files_extract \
        --model-executable \$model_executable \
        --n-cores $task.cpus \
        --output-dir . \
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
