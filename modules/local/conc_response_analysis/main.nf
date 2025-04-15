process CONC_RESPONSE_ANALYSIS {
    cpus params.n_cores

    input:
    path model
    tuple val(name), val(probes), path(all_probe_file)

    output:
    tuple val(name), path("${name}_fits_${task.index}.tar.gz"), emit: all_fits_files

    script:
    def probe_files = probes.collect { "./" + it + ".pkl" }.join(" ")
    def probe_files_extract = probes.collect { "Data/" + it + ".pkl" }.join(" ")
    """
    mkdir Data
    tar -zxf $all_probe_file -C Data/ $probe_files

    mkdir Samples
    mkdir Fits
    run_conc_response_analysis.py \
        --data-files $probe_files_extract \
        --model-name $model \
        --n-cores $task.cpus

    tar -czf ${name}_fits_${task.index}.tar.gz -C Fits/ .
    """ 
} 