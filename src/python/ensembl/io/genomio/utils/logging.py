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
from typing import Optional


# default logging formats
LOGGING_FORMAT = "%(asctime)s\t%(levelname)s\t%(message)s"
DATE_FORMAT = r"%Y-%m-%d_%H:%M:%S"

# helper functions
def _prepare_console_handler(*, debug: bool = False, verbose: bool = False) -> logging.StreamHandler:
    """setup console handler with different logging levels"""
    console_h = logging.StreamHandler()
    # set default logging level
    if debug:
        console_h.setLevel(logging.DEBUG)
    elif verbose:
        console_h.setLevel(logging.INFO)
    return console_h


def _prepare_file_handler(filename: Optional[str] = None, *, debug: bool = False) -> Optional[logging.FileHandler]:
    """setup file handler with default loggin.INFO level"""
    """retuns None if no file name defined"""
    if filename is None:
        return None

    file_h = logging.FileHandler(filename)

    file_h.setLevel(logging.INFO)
    if debug:
        file_h.setLevel(logging.DEBUG)

    return file_h


def setup_logging(
        args: argparse.Namespace,
        *, # no positional arguments allowed after this
        name: Optional[str] = None,
    ) -> logging.Logger:
    """Setup logging infrustucture."""
    """args: argparse.Namespace -- args with "debug", "verbose" and "logfile" options."""
    # get logger with a specified name (or the default one)
    logger = logging.getLogger(name)

    # set up shared formatter
    formatter = logging.Formatter(fmt = LOGGING_FORMAT, datefmt = DATE_FORMAT)

    # adding console handler
    console_h = _prepare_console_handler(debug = args.debug, verbose = args.verbose)
    console_h.setFormatter(formatter)
    logger.addHandler(console_h)

    # adding logging to file
    file_h = _prepare_file_handler(args.logfile)
    if file_h:
        file_h.setFormatter(formatter)
        logger.addHandler(file_h)

    return logger
