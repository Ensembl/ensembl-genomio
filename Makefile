# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

.PHONY: docs apidoc coverage clean

apidoc:
	sphinx-apidoc -o docs/reference/ src/python/ensembl --force --module-first --no-toc --implicit-namespaces
	rm -f docs/reference/ensembl.rst

docs: apidoc
	sphinx-build -b html docs/ docs/_build/html

coverage:
	coverage run -m pytest
	coverage html -d docs/_static/htmlcov
	coverage xml -o reports/coverage.xml
	genbadge coverage -i reports/coverage.xml -o docs/_static/htmlcov/coverage-badge.svg
	rm -rf reports

clean:
	rm -rf docs/_build docs/reference docs/_static/htmlcov
