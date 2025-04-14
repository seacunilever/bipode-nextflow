process PREPARE_INPUTS {
    publishDir "${params.results_dir}/", mode: "copy"

    input:
    path meta_data
    path meta_mapper
    path counts
    path config

    output:
    path "bifrost_inputs/*.json", emit: prepared_inputs

    script:
    """
    prepare_bifrost_inputs.py \
        --meta-data $meta_data \
        --meta-mapper $meta_mapper \
        --counts $counts \
        --config $config \
        --output-dir bifrost_inputs
    """
}
 