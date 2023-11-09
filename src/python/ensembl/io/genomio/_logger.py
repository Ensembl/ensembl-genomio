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
"""TODO"""

import logging
from pathlib import Path
from typing import Optional


def init_logging(log_level: str = "WARNING", log_file: Optional[Path] = None) -> logging.Logger:
    """TODO"""
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s", datefmt=r"%Y-%m-%d_%H:%M:%S")

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(log_level)
    stream_handler.setFormatter(formatter)
    logging.root.handlers[0] = stream_handler

    if log_file:
        logfile_handler = logging.FileHandler(log_file)
        logfile_handler.setLevel(logging.DEBUG)
        logfile_handler.setFormatter(formatter)
        logging.root.addHandler(logfile_handler)
        logging.root.setLevel(logging.DEBUG)
