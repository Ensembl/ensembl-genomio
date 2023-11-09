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
"""Utils to deal with logging."""

__all__ = ["setup_logging"]

import argparse
import logging

LOGGING_FORMAT = "%(asctime)s\t%(levelname)s\t%(message)s"
DATE_FORMAT = r"%Y-%m-%d_%H:%M:%S"

def setup_logging(args: argparse.Namespace):
    """Setup logging infrustucture."""
    """args: argparse.Namespace -- args with "debug" and "verbose" options."""
    log_level = None
    if args.debug:
        log_level = logging.DEBUG
    elif args.verbose:
        log_level = logging.INFO

    # reload(logging)
    logging.basicConfig(
        format=LOGGING_FORMAT,
        datefmt=DATE_FORMAT,
        level=log_level,
    )
