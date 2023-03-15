#!env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import eHive
import shutil
import gzip


class ungzip(eHive.BaseRunnable):
    """Simple runnable to uncompress a gzip file, and bring it as a parameter for the pipeline.

    Args:
        input: Path to the gzipped file to uncompress.
        output: Path to the final uncompressed file.
        name: Key for the dataflow.

    Dataflows:
        2: One single value for the uncompressed file, with key name.

    Note:
        If the file is not gzipped, then it will only be renamed.
    """

    def run(self):
        in_file = self.param_required("input")
        out_file = self.param_required("output")
        name = self.param("out_name")

        if in_file.endswith(".gz"):
            with gzip.open(in_file, "rb") as f_in:
                with open(out_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy(in_file, out_file)

        self.dataflow({name: out_file}, 2)
