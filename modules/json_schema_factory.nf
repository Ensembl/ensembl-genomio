process JSON_SCHEMA_FACTORY {
    tag "$manifest_dir.name"
    label 'rc_default'
    // debug true

    input:
    path manifest_dir

    output:
    path "*.json", includeInputs: true

    script:
    // Add quotes around each key of the dictionary to make the list compatible with Bash
    metadata_types = "['" + params.json_schemas.keySet().join("','") + "']"
    """
    json_schema_factory --manifest ${manifest_dir} --metadata_types "${metadata_types}"
    """
}
