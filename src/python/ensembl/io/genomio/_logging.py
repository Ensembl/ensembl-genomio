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
"""Internal logging system for the Ensembl GenomIO modules."""

import logging
from pathlib import Path
from typing import Optional, Union


def init_logging(log_level: Union[int, str] = "WARNING", log_file: Optional[Path] = None) -> None:
    """Initialises the logging system.

    By default, all the log messages corresponding to `log_level` (and) above will be printed in the standard
    error. If `log_file` is provided, all debug messages (and above) will be written into the file.

    Args:
        log_level: Minimum logging level for the standard error.
        log_file: Logging file where to write debug (and above) logging messages.

    """
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s", datefmt=r"%Y-%m-%d_%H:%M:%S")
    # Set the standard error logging handler and replace the existing one in the root logging system
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(log_level)
    stream_handler.setFormatter(formatter)
    logging.root.handlers[0] = stream_handler
    if log_file:
        # Set the debug-level logging file handler and add it to the root logging system
        logfile_handler = logging.FileHandler(log_file)
        logfile_handler.setLevel(logging.DEBUG)
        logfile_handler.setFormatter(formatter)
        logging.root.addHandler(logfile_handler)
        logging.root.setLevel(logging.DEBUG)
