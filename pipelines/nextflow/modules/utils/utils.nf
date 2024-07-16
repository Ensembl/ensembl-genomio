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

def generate_url(protocol, host, port, user, password=null, database=null) {
    // Generate a URL when all the credentials are provided
    base_url = "${protocol}://${user}@${host}:${port}"

    if (password) {
        base_url = "${protocol}://${user}:${password}@${host}:${port}"
    }
    if (database) {
        base_url += "/${database}"
    }
    return base_url
}

def extract_url_args(url_string) {
    // Extract the MySQL arguments from the URL
    def pattern = ~/mysql:\/\/(.*?)(:(.*?)?)?@(.*?):(\d+)\/?(.*)?/

    if (!(url_string =~ pattern)) {
      return [error: "Invalid url ${url_string}. A simple Mysql URL is expected."]
    }
    def url_parts = url_string.tokenize('/@')
    def user_info = url_parts[1].split(":")
    def user = user_info[0]
    def password = user_info.length > 1 ? user_info[1] : null
    def (host, port) = url_parts[2].split(":")
    def database = url_parts.size() > 3 ? url_parts[3] : null

    return [user: user, password: password, host: host, port: port, database: database]
}
