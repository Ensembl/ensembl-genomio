// See the NOTICE file distributed with this work for additional information
//regarding copyright ownership.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
"""Utils to deal with nextflow parameters."""

def Param_multi_values(String multi_value = "", List<String> allowed_values, Boolean run_all = false) {
    """
    Support parsing of multi-value parameters.
    """
    if (run_all) {
        return allowed_values
    }
    if (multi_value.isEmpty()){
        return "No value specified $multi_value"
    }
    else {
        all_params = multi_value.tokenize(',')
        all_params.each { value ->
            if (!allowed_values.contains(value)) {
                throw new Exception("Invalid param value: $value.
                    Allowed values are: $allowedValues")
            }
        }
    }
    return all_params
}

