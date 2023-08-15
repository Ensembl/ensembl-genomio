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

process EXTRACT_GENE_LISTS {
    label 'local'

    input:
        path(events, stageAs: "events.tab")

    output:
        path("new_genes.tab")
        path("changed_genes.tab")

    script:
    def changed_genes = "changed_genes.tab"
    def new_genes = "new_genes.tab"
    """
    grep -v "//" $events | grep -v "=" | grep -v "~" | sed s'/[+><]/\\t/' | cut -f3 | sed 's/:/\\n/g' | sort -u > $new_genes
    grep "=" $events | cut -f2 | sed -r 's/=[+!-]?/\t/' > $changed_genes
    """
}
