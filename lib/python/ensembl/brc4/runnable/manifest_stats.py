#!env python3

import eHive
import gzip
import sys, io, re, os
import json
from statistics import mean
from BCBio.GFF import GFFExaminer


class manifest_stats(eHive.BaseRunnable):

    def run(self):
        manifest_path = self.param_required("manifest")
        manifest = self.get_manifest(manifest_path)
        
        stats = []
        if "gff3" in manifest:
            stats += self.get_gff3_stats(manifest["gff3"])
        if "seq_region" in manifest:
            stats += self.get_seq_region_stats(manifest["seq_region"])
        
        stats_path = os.path.join(os.path.dirname(manifest_path), "stats.txt")
        with open(stats_path, "w") as stats_out:
            stats_out.write("\n".join(stats))
        
    def get_manifest(self, manifest_path):

        with open(manifest_path) as manifest_file:
            manifest = json.load(manifest_file)
            
            # Use dir name from the manifest
            for name in manifest:
                if "file" in manifest[name]:
                    file_name = manifest[name]["file"]
                    file_name = os.path.join(os.path.dirname(manifest_path), file_name)
                    manifest[name] = file_name
                else:
                    for f in manifest[name]:
                        if "file" in manifest[name][f]:
                            file_name = manifest[name][f]["file"]
                            file_name = os.path.join(os.path.dirname(manifest_path), file_name)
                            manifest[name][f] = file_name
            
            return manifest
    def get_seq_region_stats(self, seq_region_path):
        
        seq_regions = self.get_json(seq_region_path)
        
        # Get basic data
        coord_systems = {}
        circular = 0
        locations = []
        for seqr in seq_regions:
            coord_level = seqr["coord_system_level"]
            if not coord_level in coord_systems:
                coord_systems[coord_level] = []
            coord_systems[coord_level].append(seqr["length"])
            
            if "circular" in seqr:
                circular += 1
            if "location" in seqr:
                locations.append(seqr["location"])
        
        # Stats
        stats = []
        stats.append(os.path.basename(seq_region_path))
        stats.append("Total coord_systems %d" % len(coord_systems))
        for coord_name, lengths in coord_systems.items():
            stats.append("\nCoord_system: %s" % coord_name)
            
            stat_counts = {
                    "Total sequences" : len(lengths),
                    "Total sequence length" : sum(lengths),
                    "Minimum sequence length" : min(lengths),
                    "Maximum sequence length" : max(lengths),
                    "Mean sequence length" : mean(lengths),
                    }
            for name, count in stat_counts.items():
                stats.append("%9d\t%s" % (count, name))
        
        # Special
        if circular or locations:
            stats.append("\nSpecial")
            if circular:
                stats.append("%9d\t%s" % (circular, "circular sequences"))
            if locations:
                stats.append("%9d\t%s" % (len(locations), "sequences with location"))
                for loc in locations:
                    stats.append("\t\t\t%s" % loc)
        
        stats.append("\n")

        return stats

    def get_json(self, json_path):
        with open(json_path) as json_file:
            return json.load(json_file)

    def get_gff3_stats(self, gff3_path):
        
        stats = []
        stats.append(os.path.basename(gff3_path))
        if gff3_path.endswith(".gz"):
            with io.TextIOWrapper(gzip.open(gff3_path, "r")) as gff3_handle:
                stats += self.parse_gff3(gff3_handle)
        else:
            with open(gff3_path, "r") as gff3_handle:
                stats += self.parse_gff3(gff3_handle)
        stats.append("\n")
        
        return stats

    def parse_gff3(self, gff3_handle):
        biotypes = {}

        examiner = GFFExaminer()
        inside = examiner.available_limits(gff3_handle)
        print()
        
        gff_data = inside["gff_source_type"]
        
        stats = []
        for source in gff_data.keys():
            count = gff_data[source]
            biotype = source[1]
            stats.append("%9d\t%s" % (count, biotype))
        
        return stats

