"""Support module for pubmed > bibtex converter.

Contacts Entrez via BioPython efetch and esearch functions and returns
dictionary objects of the information needed for a bibtex format.

Highly dependent on python dictionary functionality, which has a tendancy
to globalise seemingly local variables as dictionary assignment creates links
between objects rather than creating new objects.
"""

from Bio import Entrez


def esearch(term, retmax):
    """Search pubmed database via Entrez.

    Expects string for term and int for retmax, as used by esearch
    Returns a list of ids of papers matching search term,
    see Biopython Entrez documentation for more details.
    """
    handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()
    return result["IdList"]


def get_record(number="18206625"):
    """Fetch and parse XML records from pubmed.

    Given a pubmed ID number as a string, return a list of the parsed xml
    - requires internet access.
    """
    handle = Entrez.efetch("pubmed", id=number, retmode="xml")
    records = Entrez.parse(handle)
    a = list(records)
    handle.close()
    return a


def summarise_record(records):
    """Short-cut function.

    Prints the title of any papers returned by get_record()
    useful for development to check that the correct results have been
    returned.
    """
    for record in records:
        print(record['MedlineCitation']['Article']['ArticleTitle'])


def list_expand(records, st, dct):
    """List version of expand function.

    Recursively calls itself and dict version to fully investigate
    tree of list and dict items (records), passing a dictionary
    which is populated in reference to the keys of the reference dct input.
    """
    if hasattr(records, '__iter__'):
        if isinstance(records, list):
            for item in records:
                dct = list_expand(item, st, dct)
            return dct
        else:
            for item in records:
                dct = dict_expand(item, records[item], st, dct)
            return dct
    else:
        return dct


def dict_expand(key, records, st, dct):
    """Dict version of expand function.

    Recursively calls itself and the list version to fully
    investigate tree of list and dict items(records), passing a dictionary
    which is populated in reference to the keys of the reference dct input.

    st is a string variable which tracks the dict keys used to get to this
    position in the tree structure.

    Printing st shows a tab separated table showing the path to variables.
    """
    st = st + key + "\t"
    tdct = check(key, records, dct)
    if hasattr(records, '__iter__'):
        if isinstance(records, list):
            for item in records:
                tdct = list_expand(item, st, tdct)
            return tdct
        else:
            for item in records:
                tdct = dict_expand(item, records[item], st, tdct)
            return tdct
    else:
        return tdct


def check(key, entry, rdict):
    """Check if a dictionary item is in a list of keywords.

    Called by dict_sort() to test if this entry needs to be added to the
    output dictionary.
    """
    if key in rdict.keys():
        rdict[key] = entry
    return rdict


def from_entrez(reference_dictionary,
                entrez_term, email, retmax=20, search=False):
    """Main called function.

    Gets records from entrez as defined by ids or search strings
    and returns a list of dictionaries with information required
    by the converter module.

    reference_dictionary expects the global variable ref_dict which defines
    the terms that will be selected from the parsed XML

    entrez_term is a comma seperated string of PMIDs ("26944449,18206625")
    or a search string, to be used with search=False and True.
    """
    set_email(email)            # sets global Entrez.email variable
    if search:
        entrez_list = esearch(entrez_term, retmax)
        entrez_term = ",".join(entrez_list)  # get and format pmids
    records = get_record(entrez_term)        # corresponding to search string
    output = []
    for record in records:
        output.append(list_expand(record, "", reference_dictionary.copy()))
    return output

ref_dict = {"AuthorList": "",
            "PubDate": "",
            "Volume": "",
            "PublicationTypeList": "",
            "ArticleTitle": "",
            "Title": "",        # Journal title
            "Issue": "",        # number
            "MedlinePgn": "",
            "PMID": "",
            "PublisherName": "",
            "AbstractText": ""
            }

def set_email(email):
    """Always tell Entrez who you are"""
    Entrez.email = email

Entrez.tool = "pubmed-to-bibtex converter"

if __name__ == "__main__":
    r = from_entrez(ref_dict, "heatmap", retmax=3, search=True)
    print r
