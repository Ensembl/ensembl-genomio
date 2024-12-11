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
"""Creates a manifest file in a folder depending on the file names ends."""

import ensembl.io.genomio
from ensembl.io.genomio.manifest.manifest import Manifest
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def main() -> None:
    """Main entrypoint."""
    parser = ArgumentParser(
        description="Compare the genomic data between the files present in a manifest file."
    )
    parser.add_argument_dst_path(
        "--manifest_dir", required=True, help="Folder where to create a manifest file"
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments()
    args = parser.parse_args()
    init_logging_with_args(args)

    manifest = Manifest(args.manifest_dir)
    manifest.create()


if __name__ == "__main__":
    main()
