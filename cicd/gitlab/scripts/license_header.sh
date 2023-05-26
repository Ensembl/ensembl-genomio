#!/usr/bin/env sh
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

FILE_EXT="\(pl\|pm\|py\|nf\|config\|\(my\|pg\|\)sql\|sqlite\|bash\|sh\|toml\|yml\)"

LICENSE_HEADER=(
    'See the NOTICE file distributed with this work for additional information'
    'regarding copyright ownership\.'
    ''
    'Licensed under the Apache License, Version 2\.0 \(the "License"\);'
    'you may not use this file except in compliance with the License\.'
    'You may obtain a copy of the License at'
    ''
    'http://www\.apache\.org/licenses/LICENSE-2\.0'
    ''
    'Unless required by applicable law or agreed to in writing, software'
    'distributed under the License is distributed on an "AS IS" BASIS,'
    'WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied\.'
    'See the License for the specific language governing permissions and'
    'limitations under the License\.'
)

PATTERN=""
for line in "${LICENSE_HEADER[@]}"; do
    if [ -z "$line" ]; then
        PATTERN="${PATTERN}[#-/]*\\s*"
    else
        PATTERN="${PATTERN}[#-/]*\\s*${line}\\n"
    fi
done

find . -regex ".*\.$FILE_EXT" -print0 | xargs -0 -I {} grep -LPzo "$PATTERN" {}
