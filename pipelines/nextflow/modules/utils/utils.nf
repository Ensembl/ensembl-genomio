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
