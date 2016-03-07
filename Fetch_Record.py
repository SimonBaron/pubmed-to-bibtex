"""Support module for pubmed > bibtex converter.

Contacts Entrez via BioPython efetch function and returns a
dictionary object of the information needed for a bibtex format
"""

from Bio import Entrez

def get_record(number):
    """Given a pubmed ID number, return a list of the parsed xml
    - requires internet access and closes connection"""
    handle = Entrez.efetch("pubmed", id=number, retmode="xml")
    records = Entrez.parse(handle)
    a = list(records)
    handle.close()
    return a


def summarise_record(records):
    """Short-cut function prints the title of any papers returned
    by get_record()"""
    for record in records:
        print(record['MedlineCitation']['Article']['ArticleTitle'])

def list_expand(records, st, dct):
    """List version of expand function.

    Recursively calls itself and dict version to fully investigate
    tree of list and dict items, passing a dictionary
    which is populated in reference to the keys of the referenct dict input"""
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
    """Dict version of expand function

    Recursively calls itself and the list version to fully
    investigate tree of list and dict items, passing a dictionary
    which is populated in reference to the keys of the referenct dict input
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


def from_entrez(reference_dictionary, entrez_id):
    """Main called function.

    Gets records from entrez as defined by ids and returns a list of
    dictionaries with information requred by the converter module.
    """
    records = get_record(entrez_id)
    output = []
    for record in records:
        output.append(list_expand(record, "", reference_dictionary))
    return output

ref_dict = {"AuthorList": "",
            "DateCompleted": "",
            "Volume": "",
            "PublicationTypeList": "",
            "ArticleTitle": "",
            "Title": "",        # Journal title
            "Issue": "",        # number
            "MedlinePgn": ""
            }

Entrez.email = "simon.c.baron@gmail.com"
