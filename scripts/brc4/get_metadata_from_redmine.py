#!env python3

from redminelib import Redmine
import argparse
import os, json
 
url = 'https://redmine.apidb.org'
default_fields = dict(
        status_id = 'open',
        cf_17 = "Data Processing (EBI)",
        )

def retrieve_genomes(redmine, output_dir, build=None):
    """
    Get genomes metadata from Redmine, store them in json files.
    Each issue/new_genome is stored as one file in the output dir
    """
    
    issues = get_genome_issues(redmine, build)
    if not issues:
        print("No issues found")
        return
    else:
        print("%d issues found" % len(issues))
    
    # Create the output dir
    try:
        os.mkdir(output_dir)
    
    for issue in issues:
        genome_structure = parse_genome(issue)
        organism = str(issue.id)

        organism_file = output_dir + "/" + organism + ".json"
        f = open(organism_file, "w")
        json.dump(genome_structure, f, indent=True)
        f.close()
        
def parse_genome(issue):
    """
    Extract genome metadata from a Redmine issue
    Return a nested dict
    """
    
    customs = get_custom_fields(issue)
    genome = {
            "BRC4": {},
            "species": {},
            "assembly": {},
            "genebuild": {},
            }
    
    if "GCA number" in customs:
        genome["assembly"]["accession"] = customs["GCA number"]["value"]
    if "Component DB" in customs:
        components = customs["Component DB"]["value"]
        if len(components) == 1:
            genome["BRC4"]["component"] = components[0]
        elif len(components) > 1:
            raise Exception("More than 1 component for new genome " + issue.name)

    return genome

def get_custom_fields(issue):
    """
    Put all Redmine custom fields in a dict instead of an array
    Return a dict
    """
    
    cfs = {}
    for c in issue.custom_fields:
        cfs[c["name"]] = c
    return cfs

def get_genome_issues(redmine, build=None):
    """Retrieve all issue for new genomes, be they with or without gene sets"""
    
    other_fields = dict(
            cf_94 = [
                "Genome sequence and Annotation",
                "Assembled genome sequence without annotation"
                ],
            )
    if build:
        other_fields["fixed_version"] = 'Build ' + str(build)

    return get_ebi_issues(redmine, other_fields)
    
def get_ebi_issues(redmine, other_fields=dict()):
    """
    Get EBI issues from Redmine, add other fields if provided
    Return a Redmine ResourceSet
    """
    
    search_fields = { **default_fields, **other_fields }
    
    return redmine.issue.filter(**search_fields)


def main():
    parser = argparse.ArgumentParser(description='Retrieve metadata from Redmine')
    
    parser.add_argument('--key', type=str, required=True,
                help='Redmine authentification key')
    parser.add_argument('--output_dir', type=str, required=True,
                help='Output_dir')
    # Choice
    parser.add_argument('--get', choices=['genomes', 'rnaseq', 'dnaseq'], required=True,
                help='Get genomes, rnaseq, or dnaseq issues')
    # Optional
    parser.add_argument('--build', type=int,
                help='Restrict to a given build')

    args = parser.parse_args()
    
    redmine = Redmine(url, key=args.key)
    if args.get == 'genomes':
        retrieve_genomes(redmine, args.output_dir, args.build)
    elif args.get == 'rnaseq':
        print("RNA-Seq Redmine retrieval to be implemented")
    elif args.get == 'dnaseq':
        print("DNA-Seq Redmine retrieval to be implemented")
    else:
        print("Need to say what data you want to get: --genomes? --rnaseq? --dnaseq?")

if __name__ == "__main__":
    main()
