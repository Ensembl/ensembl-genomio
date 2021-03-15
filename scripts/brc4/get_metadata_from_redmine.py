#!env python3

from redminelib import Redmine
import argparse
import os, json, re
 
url = 'https://redmine.apidb.org'
default_fields = dict(
        status_id = 'open',
        cf_17 = "Data Processing (EBI)",
        )
insdc_pattern = "^GC[AF]_\d{9}(\.\d+)?$"

def retrieve_genomes(redmine, output_dir, build=None):
    """
    Get genomes metadata from Redmine, store them in json files.
    Each issue/new_genome is stored as one file in the output dir
    """
    
    genomes_with_genes = get_issues(redmine, "Genome sequence and Annotation", build)
    genomes_without_genes = get_issues(redmine, "Assembled genome sequence without annotation", build)
    
    issues = genomes_with_genes + genomes_without_genes
    if not issues:
        print("No issues found")
        return
    else:
        print("%d issues found" % len(issues))
    
    # Create the output dir
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    for issue in issues:
        genome_structure = parse_genome(issue)
        if not genome_structure: continue
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
    
    # Get GCA accession
    if "GCA number" in customs:
        accession = customs["GCA number"]["value"]
        accession = check_accession(accession)
        if not accession:
            print("No proper accession for issue %d (%s): %s" % (issue.id, issue.subject, customs["GCA number"]["value"]))
            return
        genome["assembly"]["accession"] = accession
    else:
        print("No accession for issue %d (%s)" % (issue.id, issue.subject))

    # Get BRC4 component
    if "Component DB" in customs:
        components = customs["Component DB"]["value"]
        if len(components) == 1:
            genome["BRC4"]["component"] = components[0]
        elif len(components) > 1:
            raise Exception("More than 1 component for new genome " + issue.name)
    else:
        print("No component for issue %d (%s)" % (issue.id, issue.subject))

    return genome

def check_accession(accession):
    """
    Check the accession string format
    Returns a cleaned up accession
    """
    accession = accession.strip()
    
    # Remove the url if it's in one
    # There might even be a trailing url
    accession = re.sub(r"^.+/([^/]+)/?", r"\1", accession)
    
    if re.match(insdc_pattern, accession):
        return accession
    else:
        return

def get_custom_fields(issue):
    """
    Put all Redmine custom fields in a dict instead of an array
    Return a dict
    """
    
    cfs = {}
    for c in issue.custom_fields:
        cfs[c["name"]] = c
    return cfs

def get_issues(redmine, datatype, build=None):
    """
    Retrieve all issue for new genomes, be they with or without gene sets
    Return a Redmine ResourceSet
    """
    
    other_fields = { "cf_94" : datatype }
    if build:
        other_fields["fixed_version"] = 'Build ' + str(build)

    return list(get_ebi_issues(redmine, other_fields))
    
def get_ebi_issues(redmine, other_fields=dict()):
    """
    Get EBI issues from Redmine, add other fields if provided
    Return a Redmine ResourceSet
    """

    # Other fields replace the keys that already exist in default_fields
    search_fields = { **default_fields, **other_fields }
    
    return redmine.issue.filter(**search_fields)


def main():
    # Parse command line arguments
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
    
    # Start Redmine API
    redmine = Redmine(url, key=args.key)
    
    # Choose which data to retrieve
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
