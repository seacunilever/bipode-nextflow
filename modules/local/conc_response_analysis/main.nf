process CONC_RESPONSE_ANALYSIS {
    tag "${meta.id}"
    cpus params.n_cores

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

    input:
    path model
    tuple val(meta), val(probes), path(all_probe_file)

    output:
    tuple val(meta), path("${prefix}.tar.gz"), emit: compressed_fits_files

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def probe_files = probes.collect { "./" + it + ".pkl" }.join(" ")
    def probe_files_extract = probes.collect { "Data/" + it + ".pkl" }.join(" ")
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir Data
    tar -zxf $all_probe_file -C Data/ $probe_files

    mkdir Samples
    mkdir Fits
    run_conc_response_analysis.py \
        --data-files $probe_files_extract \
        --model-name $model \
        --n-cores $task.cpus \
        $args

    sleep 5

    tar $args2 -cf - -C Fits/ . | gzip $args3 > ${prefix}.tar.gz
    """
}
