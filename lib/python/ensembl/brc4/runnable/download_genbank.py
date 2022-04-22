#!env python3

import eHive
import os, re
import requests

class download_genbank(eHive.BaseRunnable):
    
    def param_defaults(self):
        return {
        }

    def run(self):
        accession = self.param_required('gb_accession')
        main_download_dir = self.param('download_dir')
        download_dir = main_download_dir + "/" + accession

        # Set and create dedicated dir for download
        if not os.path.isdir(download_dir):
            os.makedirs(download_dir)

        # Download file
        gb_path = self.download_genbank(accession, download_dir)
        
        output = { 'gb_file' : gb_path, 'gb_accession' : accession }
        self.dataflow(output, 2)


    def download_genbank(self, accession, dl_dir):
        """
        Given a GenBank accession, download the corresponding file in GenBank format
        """
        dl_path = os.path.join(dl_dir, accession + '.gb')
        if os.path.exists(dl_path):
            return dl_path

        # Get the list of assemblies for this accession
        e_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        e_params = {
                "db" : "nuccore",
                "rettype" : "gb",
                "retmode" : "text",
                }
        e_params["id"] = accession
        
        result = requests.get(e_url, params=e_params)
        
        if result and result.status_code == 200:
            with open(dl_path, "wb") as gff:
                gff.write(result.content)
            print(f"GFF file write to {dl_path}")
            return dl_path
        else:
            print("Error: " + str(result))
            
