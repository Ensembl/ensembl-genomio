process CHECK_JSON_SCHEMA {
    tag "$json_file.name"
    label 'rc_default'
    // debug true

    input:
    path json_file

    script:
    schema = params.json_schemas[json_file.baseName]
    """
    check_json_schema --json_file ${json_file} --json_schema ${schema}
    """
}