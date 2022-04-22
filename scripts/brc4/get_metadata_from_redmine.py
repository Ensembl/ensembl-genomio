#!env python3

from redminelib import Redmine
import argparse
import os, json, re, time
import requests
import xml.etree.ElementTree as ET
 
url = 'https://redmine.apidb.org'
default_fields = dict(
        status_name = 'Data Processing (EBI)',
        cf_17 = "Data Processing (EBI)",
        )
insdc_pattern = "^GC[AF]_\d{9}(\.\d+)?$"
accession_api_url = "https://www.ebi.ac.uk/ena/browser/api/xml/%s"
veupathdb_id = 1976

def retrieve_genomes(redmine, output_dir, build=None):
    """
    Get genomes metadata from Redmine, store them in json files.
    Each issue/new_genome is stored as one file in the output dir
    """
    issues = get_all_genomes(redmine, build)
    if not issues:
        print("No files to create")
        return
    
    # Create the output dir
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    ok_new = []
    ok_patch = []
    ok_other = []
    failed_issues = []
    
    groups = {
            "new_genomes" : [],
            "copy_ensembl": [],
            "other": [],
            }

    for issue in issues:
        genome, extra = parse_genome(issue)
        failure = check_genome(genome, extra)
        
        if failure:
            failed_issues.append({"issue" : issue, "desc" : failure})
            continue
    
        abbrev = genome["BRC4"]["organism_abbrev"]
        group = "other"
        if "Load from INSDC" in extra["operations"] or "Load from RefSeq" in extra["operations"]:
            group = "new_genomes"
            ok_new.append({"issue" : issue, "desc" : abbrev})
        elif "Load from EnSEMBL" in extra["operations"]:
            group = "copy_ensembl"
            ok_new.append({"issue" : issue, "desc" : abbrev})
        elif "Patch build" in extra["operations"]:
            group = "patch_build"
            ok_patch.append({"issue" : issue, "desc" : abbrev})
        else:
            group = "other"
            ok_other.append({"issue" : issue, "desc" : abbrev})
        
        groups[group].append(genome)
    
    # Write files
    for group, genomes in groups.items():
        if genomes:
            group_dir = os.path.join(output_dir, group)
            try:
                os.makedirs(group_dir)
            except FileExistsError:
                pass
                
            
            for genome in genomes:
                organism = genome["BRC4"]["organism_abbrev"]
                organism_file = os.path.join(group_dir, organism + ".json")
                with open(organism_file, "w") as f:
                    json.dump(genome, f, indent=True)
            

    # Print summaries
    print_summary(failed_issues, "failed issues")
    print_summary(ok_other, "other genome operations")
    print_summary(ok_patch, "patch builds")
    print_summary(ok_new, "new genomes")

def get_all_genomes(redmine, build=None):
    """
    Query Redmine to get all new genomes, with or without genes
    """
    genomes_with_genes = get_issues(redmine, "Genome sequence and Annotation", build)
    genomes_without_genes = get_issues(redmine, "Assembled genome sequence without annotation", build)
#    genomes_without_genes = []
    print("%d issues for new genomes with genes found" % len(genomes_with_genes))
    print("%d issues for new genomes without genes found" % len(genomes_without_genes))
    
    issues = genomes_with_genes + genomes_without_genes
    
    print("%d issues found" % len(issues))
    return issues
 
def parse_genome(issue):
    """
    Extract genome metadata from a Redmine issue
    Return a nested dict
    """
    
    customs = get_custom_fields(issue)
    genome = {
            "BRC4": {
                "component" : "",
                "organism_abbrev" : "",
                },
            "species": {},
            "assembly": {
                "accession" : ""
                },
            "genebuild": {},
            }
    
    extra = {}
    
    # Get GCA accession
    if "GCA number" in customs:
        accession = customs["GCA number"]["value"]
        accession = check_accession(accession)
        genome["assembly"]["accession"] = accession

    # Get BRC4 component
    if "Component DB" in customs:
        components = customs["Component DB"]["value"]
        if len(components) == 1:
            genome["BRC4"]["component"] = components[0]
        elif len(components) > 1:
            raise Exception("More than 1 component for new genome " + str(issue.id))

    # Get Organism abbrev
    try:
        abbrev = customs["Organism Abbreviation"]["value"]
        if abbrev:
            genome["BRC4"]["organism_abbrev"] = abbrev
    except:
        print("Can't get organism abbrev for %s" % issue.id)

    # Warn to get GFF2Load
    try:
        gff_path = customs["GFF 2 Load"]["value"]
        if gff_path:
            extra["GFF"] = True
            #print("GFF2LOAD: separate gff file for %s: %s (issue %d)" % (genome["BRC4"]["organism_abbrev"], gff_path, issue.id))
    except:
        pass

    # Warn for replacement
    try:
        if customs["Replacement genome?"]["value"].startswith("Yes"):
            extra["Replacement"] = True
            #print("REPLACEMENT: the organism %s is a replacement (issue %d)" % (genome["BRC4"]["organism_abbrev"], issue.id))
    except:
        pass

    # Operations
    try:
        extra["operations"] = get_operations(issue)
        
    except:
        pass

    return (genome, extra)

def check_genome(genome, extra):
    
    if not genome:
        return "No genome parsed"
    
    if not "organism_abbrev" in genome["BRC4"] or not genome["BRC4"]["organism_abbrev"]:
        return "No organism_abbrev defined"
    
    operations = extra["operations"]
    
    if not operations:
        return "No EBI operation defined"
    else:
        source_count = 0
        if "Load from INSDC" in operations: source_count += 1
        if "Load from RefSeq" in operations: source_count += 1
        if "Load from EnsEMBL" in operations: source_count += 1
        
        if source_count > 1:
            return f"Mix of {source_count} sources"
        if source_count == 1:
            if not "accession" in genome["assembly"] or not genome["assembly"]["accession"]:
                return "No accession for assembly"
        else:
            # Other operations
            if "Other" in "operations":
                return ""
            
    
    return ""

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
        version_id = get_version_id(redmine, build)
        other_fields["fixed_version_id"] = version_id

    return list(get_ebi_issues(redmine, other_fields))

def get_version_id(redmine, build):
    """
    Given a build number, get the version id for it
    """
    versions = redmine.version.filter(project_id=veupathdb_id)
    version_name = "Build " + str(build)
    version_id = [version.id for version in versions if version.name == version_name]
    return version_id
    
def get_ebi_issues(redmine, other_fields=dict()):
    """
    Get EBI issues from Redmine, add other fields if provided
    Return a Redmine ResourceSet
    """

    # Other fields replace the keys that already exist in default_fields
    search_fields = { **default_fields, **other_fields }
    
    return redmine.issue.filter(**search_fields)

def add_genome_organism_abbrev(redmine, build, abbrevs_file, update=False):
    """
    Retrieve genome issues, get the corresponding INSDC accession, and generate an organism abbrev
    Put the abbrev back to the Redmine ticket
    By default do not update the ticket
    """
    
    all_abbrevs = load_abbrevs(abbrevs_file)
    
    # First get the genomes issues
    issues = get_all_genomes(redmine, build)
    if not issues:
        print("No Redmine tickets to update")
        return
    
    # Keep track of each problem
    ok_generate = []
    ok_exist = []
    ok_replace = []
    warnings_replace = []
    warnings = []
    failed_issues = []
    not_new_ok = []
    not_new_fail = []
    
    for issue in issues:
        time.sleep(0.1)
        (genome, extra) = parse_genome(issue)
        custom = get_custom_fields(issue)
        if not genome:
            failed_issues.append({"issue": issue, "desc": "Not enough information to parse"})
            continue
        
        cur_abbrev = custom["Organism Abbreviation"]["value"]

        if not is_new_genome(issue):
            if cur_abbrev:
                not_new_ok.append({"issue": issue, "desc": f"Not a new genome ({cur_abbrev})"})
            else:
                not_new_fail.append({"issue": issue, "desc": "Please add organism_abbrev manually"})
            continue

        try:
            accession = genome["assembly"]["accession"]
            full_name = custom["Experimental Organisms"]["value"]
            new_abbrev = make_organism_abbrev(full_name)

            if new_abbrev in all_abbrevs:
                if "Replacement" in extra:
                    ok_replace.append({"issue" : issue, "desc" : new_abbrev})
                else:
                    failed_issues.append({
                        "issue" : issue,
                        "desc" : f"Abbrev {new_abbrev} used by other species"})
            elif cur_abbrev:
                if cur_abbrev == new_abbrev:
                    if "Replacement" in extra and not new_abbrev in all_abbrevs:
                        warnings_replace.append({"issue" : issue, "desc" : f"replacement genome has new abbrev: {new_abbrev}"})
                    else:    
                        ok_exist.append({"issue" : issue, "desc" : new_abbrev})
                else:
                    warnings.append({
                        "issue" : issue,
                        "desc" : f"Generated id differs: {cur_abbrev} vs {new_abbrev}"})
                    
            else:
                if "Replacement" in extra and not new_abbrev in all_abbrevs:
                    warnings_replace.append({"issue" : issue, "desc" : f"replacement genome has new abbrev: {new_abbrev}"})
                else:    
                    ok_generate.append({"issue" : issue, "desc" : new_abbrev})
                add_organism_to_issue(redmine, issue, new_abbrev, update)
        except Exception as e:
            failed_issues.append({"issue" : issue, "desc" : f"Can't make: {e}"})
            continue
    
    print("\n%d issues considered" % len(issues))

    print_summary(ok_generate, "issues with generated organism abbrev (use --update to set the field)")
    print_summary(ok_exist, "issues with existing organism_abbrev")
    print_summary(ok_replace, "issues are replacement with existing organism_abbrev")
    print_summary(warnings_replace, "issues are replacement but without a previous organism_abbrev")
    print_summary(warnings, "issues with warnings")
    print_summary(failed_issues, "issues failed")
    print_summary(not_new_ok, "not new genomes, with organism_abbrev")
    print_summary(not_new_fail, "not new genomes, need organism_abbrev")

def print_summary(summaries, description):
    if summaries:
        print()
        print(f"{len(summaries)} {description}:")
        for summary in summaries:
            desc = summary["desc"]
            issue = summary["issue"]
            operations = get_operations(issue)
            gff = " +GFF" if has_gff(issue) else ""
            stable_ids = " +STABLE_IDS" if has_stable_ids(issue) else ""
            replace = " +REPLACE" if is_replacement(issue) else ""
            ops = ",".join(operations)
            desc = f"{desc} ({ops}{gff}{stable_ids}{replace})"
            print(f"\t{desc:64}\t{issue.id:8}\t{issue.subject}")

def get_operations(issue):
    customs = get_custom_fields(issue)
    return customs["EBI operations"]["value"]

def is_new_genome(issue):
    operations = get_operations(issue)
    if "Load from INSDC" in operations or "Load from RefSeq" in operations or "Load from EnsEMBL" in operations:
        return True
    else:
        return False

def is_patch_build(issue):
    customs = get_custom_fields(issue)
    if customs["Patch build"]["value"]:
        return True
    else:
        return False

def is_replacement(issue):
    customs = get_custom_fields(issue)
    if customs["Replacement genome?"]["value"].startswith("Yes"):
        return True
    else:
        return False

def has_gff(issue):
    customs = get_custom_fields(issue)
    if customs["GFF 2 Load"]["value"]:
        return True
    else:
        return False

def has_stable_ids(issue):
    operations = get_operations(issue)
    if "Allocate stable ids" is operations:
        return True
    else:
        return False

def add_organism_to_issue(redmine, issue, new_abbrev, update=False):
    """
    Actually update the Redmine tickets with the organism_abbrev
    """
    if update:
        print("Add organism Abbrev %s for issue %d" % (new_abbrev, issue.id))
        custom = get_custom_fields(issue)
        
        if "Organism Abbreviation" in custom:
            field_id = custom["Organism Abbreviation"]["id"]
            feedback = redmine.issue.update(issue.id,
                    custom_fields=[{ 'id' : field_id, 'value' : new_abbrev }]
                    )
            if not feedback:
                print("Failed to update the organism abbrev (id %d) with %s" % (field_id, new_abbrev))
        else:
            raise Exception("Can't find organism abbrev custom key id")
    else:
        print("[DRY RUN] Add organism abbrev %s for issue %d" % (new_abbrev, issue.id))

def load_abbrevs(path):
    """
    Load a list of organism abbrevs from a file. Expected to be one per line
    """
    if not path:
        print("Warning: I don't have a list of older abbrevs to compare with.")
        return
    abbrevs = []
    with open(path, "r") as abbr_file:
        for line in abbr_file:
            line = line.rstrip()
            if line:
                fields = line.split("\t")
                if len(fields) == 1:
                    abbrevs.append(line)
                else:
                    raise Exception("Can't load current abbrevs from a multicolumn string")
    return abbrevs
    
def make_organism_abbrev(name):
    
    name = name.strip()
    if name == "":
        raise Exception("field 'Experimental Organisms' needed")
    items = name.split(" ")
    if len(items) < 3:
        raise Exception(f"name is too short ({name})")

    genus = items[0]
    species = items[1]
    strain_abbrev = "".join(items[2:])
    
    genus = re.sub("[\[\]]", "", genus)
    strain_abbrev = re.sub(r"(isolate|strain|breed|str\.|subspecies|sp\.)", "", strain_abbrev, flags=re.IGNORECASE)
    strain_abbrev = re.sub(r"[\/\(\)#:]", "", strain_abbrev)
    
    organism_abbrev = genus[0].lower() + species[0:3] + strain_abbrev
    return organism_abbrev
    

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Retrieve metadata from Redmine')
    
    parser.add_argument('--key', type=str, required=True,
                help='Redmine authentification key')
    parser.add_argument('--output_dir', type=str,
                help='Output_dir')
    # Choice
    parser.add_argument('--get', choices=['genomes', 'organism_abbrev', 'abbreviate'], required=True,
                help='Get genomes, or set organism_abbrev field (use update_redmine for actually changing it)')
    # Optional
    parser.add_argument('--build', type=int,
                help='Restrict to a given build')
    parser.add_argument('--current_abbrevs', type=str,
                help='File that contains the list of current organism_abbrevs')
    parser.add_argument('--update_redmine', action='store_true', dest='update_redmine',
                help='Actually update Redmine for the organism_abbrev (dry run by default)')
    parser.add_argument('--organism', type=str,
                help='Organism name to abbreviate')
    args = parser.parse_args()
    
    # Start Redmine API
    redmine = Redmine(url, key=args.key)
    
    # Choose which data to retrieve
    if args.get == 'genomes':
        if args.output_dir:
            retrieve_genomes(redmine, args.output_dir, args.build)
        else:
            print("Need --output_dir")
            return
    elif args.get == 'organism_abbrev':
        add_genome_organism_abbrev(redmine, args.build, args.current_abbrevs, args.update_redmine)
    elif args.get == 'abbreviate':
        abbrev = make_organism_abbrev(args.organism)
        print(abbrev)
    else:
        print("Need to say what data you want to --get: genomes? organism_abbrev?")

if __name__ == "__main__":
    main()
