#!env python3

import eHive
import json
import os, re
import ftplib

class download_assembly_data(eHive.BaseRunnable):

    def run(self):
        genome_json = self.param('genome_json')
        download_dir = self.param('download_dir')

        # Get genome data
        genome_data = self.get_json(genome_json)
        accession = genome_data["assembly"]["accession"]

        # Set dedicated dir for download
        download_dir = os.path.join(download_dir, accession)
        
        if not os.path.isdir(download_dir):
            os.makedirs(download_dir)
            self.download_files(accession, download_dir)

        # TODO: md5check
        
        # Select specific files and give them a name
        files = self.get_files_selection(download_dir)

        if len(files) == 0:
            raise Exception("No file downloaded")

        files["download_dir"] = download_dir

        # Output all those files
        self.dataflow(files, 2)

    def get_json(self, json_path):
        with open(json_path) as json_file:
            json_data = json.load(json_file)
        return json_data

    def download_files(self, accession, dl_dir):
        match = re.match("GC[AF]_([0-9]{3})([0-9]{3})([0-9]{3})\.?([0-9]+)", accession)
        part1 = match.group(1)
        part2 = match.group(2)
        part3 = match.group(3)
        version = match.group(4)
        parts = (part1, part2, part3)

        # Get the list of assemblies for this accession
        f = ftplib.FTP()
        f.connect("ftp.ncbi.nlm.nih.gov")
        f.login()
        f.cwd("genomes/all/GCA/%s/%s/%s" % parts)

        files = []
        for (ftp_dir, entry) in f.mlsd():
            if re.match(accession, ftp_dir):
                f.cwd(ftp_dir)
                
                # Get all the files
                for (ftp_file, file_entry) in f.mlsd():
                    if re.match("^\.", ftp_file): continue
                    if file_entry["type"] == "dir": continue

                    local_path = os.path.join(dl_dir, ftp_file)
                    with open(local_path, 'wb') as fp:
                        f.retrbinary("RETR " + ftp_file, fp.write)

    def get_files_selection(self, dl_dir):
        files = {}
        for dl_file in os.listdir(dl_dir):
            if dl_file.endswith("assembly_report.txt"):
                files["report"] = dl_file
            elif dl_file.endswith("genomic.fna.gz"):
                files["fasta_dna"] = dl_file
            elif dl_file.endswith("genomic.gff.gz"):
                files["gff3"] = dl_file
        return files

