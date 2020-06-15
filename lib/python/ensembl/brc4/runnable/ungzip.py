#!env python3

import eHive
import json
import shutil
import gzip

class ungzip(eHive.BaseRunnable):

    def run(self):
        in_file = self.param_required("input")
        out_file = self.param_required("output")
        name = self.param("out_name")
        
        if in_file.endswith(".gz"):
            with gzip.open(in_file, 'rb') as f_in:
                with open(out_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutils.copy(in_file, out_file)
            
        self.dataflow({ name : out_file }, 2)
