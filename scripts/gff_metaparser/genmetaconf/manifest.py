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


## MANIFEST CONF ##
import gzip
import json
import os
import shutil
import subprocess as sp
import sys

from os.path import dirname, join as pj


class Manifest:
    def __init__(self, mapping, ungzip=True, always_copy=True):
        self.files = mapping
        self.ungzip = ungzip
        self.always_copy = always_copy

    def is_gz(self, name):
        return name.endswith(".gz")

    def gunzip(self, name, tag, outdir):
        if not outdir:
            print("no out_dir specified to uncompress to", file=sys.stderr)
            return None

        nogzname = name.replace(".gz", "")
        sfx = nogzname[nogzname.rfind(".") :]
        outfile = pj(outdir, tag + sfx)

        os.makedirs(outdir, exist_ok=True)
        sp.run(r"""gunzip -c %s > %s""" % (name, outfile), shell=True)
        return outfile

    def copy(self, name, tag, outdir):
        if not outdir:
            print("no out_dir specified to copy to", file=sys.stderr)
            return None

        sfx = name[name.rfind(".") :]
        outfile = pj(outdir, tag + sfx)

        os.makedirs(outdir, exist_ok=True)
        try:
            shutil.copyfile(name, outfile)
        except shutil.SameFileError:
            pass
        return outfile

    def md5sum(self, name):
        if not name:
            return
        pre = sp.check_output("md5sum %s" % name, shell=True)
        # pre = sp.check_output("md5 %s" % name, shell=True)
        (md5sum, *rest) = pre.split()
        return md5sum

    def dump(self, json_out, outdir=None):
        if outdir:
            os.makedirs(outdir, exist_ok=True)
        out = {}
        for tag, name in self.files.items():
            if not name:
                continue
            outfile = name
            if self.is_gz(name) and self.ungzip:
                print("uncompressing %s for %s" % (name, tag), file=sys.stderr)
                outfile = self.gunzip(name, tag, outdir)
            elif self.always_copy:
                print("copying %s for %s" % (name, tag), file=sys.stderr)
                outfile = self.copy(name, tag, outdir)
            print("calculating md5 for %s (%s)" % (outfile, tag), file=sys.stderr)
            md5 = self.md5sum(outfile)
            if outfile and md5:
                out[tag] = {"file": outfile, "md5sum": md5.decode()}
        if out:
            with open(json_out, "wt") as jf:
                json.dump(out, jf, indent=2)
