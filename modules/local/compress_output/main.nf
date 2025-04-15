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