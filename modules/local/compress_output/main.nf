process COMPRESS_OUTPUT {
    publishDir "${params.results_dir}/", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

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