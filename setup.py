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

"""Build script for setuptools."""
import sys

from setuptools import setup, find_namespace_packages

with open('README.md') as f:
    readme = f.read()

with open('VERSION') as f:
    version = f.read()

def import_requirements(requirements_path):
    """Import file located at the root of the repository."""
    with open(requirements_path) as file:
        return [line.rstrip() for line in file.readlines()]
setup(
    name='ensembl-genomio',
    version=version,
    packages=find_namespace_packages(where='lib/python'),
    package_dir={"": "lib/python"},
    description="Ensembl Genome IO",
    include_package_data=True,
    install_requires=import_requirements('requirements.txt'),
    long_description=readme,
    author='Ensembl',
    author_email='dev@ensembl.org',
    url='https://www.ensembl.org',
    download_url='https://github.com/Ensembl/ensembl-genomio',
    license="Apache License 2.0",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ]
)
