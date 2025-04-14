nextflow.enable.dsl=2
import groovy.json.JsonSlurper

process SPLIT_DATA {

    input:
    tuple val(name), path(input_data)

    output:
    tuple val(name), path("${name}_probes.tar.gz"), emit: all_probe_files

    script:
    """
    split_data.py --input-file $input_data --analysis-dir .
    tar -czf ${name}_probes.tar.gz -C Data/ .
    """
}

process CONC_RESPONSE_ANALYSIS {
    cpus params.n_cores

    input:
    path model
    tuple val(name), path(all_probe_file), val(probes)

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

process COMPRESS_OUTPUT {
    publishDir "${params.results_dir}/", mode: "copy"

    input:
    tuple val(name), path(all_fits_files)

    output:
    path "${name}_Results.json.zip"

    script:
    def fits_files = all_fits_files.join(" ")
    """
    mkdir Fits
    for file in $fits_files; do tar -zxf "\$file" -C Fits/; done

    compress_output.py \
        --fits-dir Fits \
        --output ${name}_Results.json.zip
    """ 
}

def get_probes(String data_file) {
    json_slurper = new JsonSlurper()
    def dataset = json_slurper.parse(file(data_file))
    def all_probes = dataset.probes
    return all_probes.collect { tuple(file(data_file).simpleName, it) }
}

def data_files = (params.data_files in List) ? params.data_files : [params.data_files]
def probes = data_files.collectMany { get_probes(it) }
def n_probes = probes.groupBy{ it[0] }.collectEntries{ k, v -> [(k):Math.ceil(v*.get(1).size()/params.n_cores.toInteger()).toInteger()] }

data_file_ch = Channel.fromList(data_files.collect { tuple(file(it).simpleName, file(it)) })
probes_ch = Channel.fromList(probes)
                   .groupTuple(size: params.n_cores.toInteger(), remainder: true, sort: true)
model_ch = Channel.fromPath(params.model_file + ".stan")

workflow {
    SPLIT_DATA(data_file_ch)
    CONC_RESPONSE_ANALYSIS(
        model_ch.first(),
        SPLIT_DATA.out.all_probe_files.combine(probes_ch, by: 0))
    COMPRESS_OUTPUT(CONC_RESPONSE_ANALYSIS.out.all_fits_files.groupTuple())
}