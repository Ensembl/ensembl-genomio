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

import hashlib

def sha256_key(name: str, repeat_class: str, repeat_type: str, seq: str | None) -> str:
    """Compute the expected SHA-256 repeat consensus key.

    Args:
        name: Repeat name.
        repeat_class: Repeat class.
        repeat_type: Repeat type.
        seq: Consensus sequence.

    Returns:
        str:Expected SHA-256 digest.

    """
    norm = "".join((seq or "").split()).upper()
    payload = f"{name}\t{repeat_class}\t{repeat_type}\t{norm}".encode()
    return hashlib.sha256(payload).hexdigest()