// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
import groovy.json.JsonSlurper

def get_key_list(dict) {
    // Add quotes around each key of the dictionary to make the list compatible with Bash
    return "['" + dict.keySet().join("','") + "']"
}

def read_json(json_path) {
    slurp = new JsonSlurper()
    json_file = file(json_path)
    text = json_file.text
    // unfortunately
    //   return slurp.parseText(text)
    // doesn't work for a single element list, we suspect lazy eval
    // symptom: instead of `[a:..., b:...]` we see the same stuff in the curly brackets `{a:..., b:...}`
    not_a_lazy_val = slurp.parseText(text)
    return not_a_lazy_val
}

def parse_list_param(String multi_value = '', List<String> allowed_values, Boolean run_all = false) {
    //Support parsing of multi-value parameters.
    if (run_all) {
        return allowed_values
    }

    if (multi_value.isEmpty()) {
        return "No value specified $multi_value"
    }
    else {
        all_params = multi_value.tokenize(',')
        all_params.each { value ->
            if (!allowed_values.contains(value)) {
                throw new Exception("Invalid param value: $value. Allowed values are: $allowed_values")
            }
        }
    }
    return all_params
}

def generate_url(host, port, user, pass=null, dbname=null) {
    base_url = "${protocol}://${user}@${host}:${port}"

    if (pass) {
        base_url = "${protocol}://${user}:${password}@${host}:${port}"
    }

    if (dbname) {
        base_url += "/${dbname}"
    }
    return base_url
}

def extractMySQLArguments(urlString) {
    def remove_protocol = urlString.split("//")[1]

    def urlparts = remove_protocol.split("@")
    def user_info = urlparts[0].split(":")
    def user = user_info[0]
    def pass = user_info.length > 1 ? user_info[1]:null

    def host_database = urlparts[1].split("/")
    def host_info = host_database[0].split(":")
    def host = host_info[0]
    def port = host_info[1]

    def database = host_database.length > 1 ? host_database[1] : null
    
    def result = [user: user, pass: pass, host: host, port: port, database: database]
    
    return result
}

