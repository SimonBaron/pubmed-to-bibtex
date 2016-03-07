'''Contact Entrez via BioPython efetch function and return a dictionary object of the information needed for a bibtex format'''

from Bio import Entrez

Entrez.email = "simon.c.baron@gmail.com"


def get_record(number):
    "Given a pubmed ID number, return a list of the parsed xml \
    - requires internet access and closes connection"

    handle = Entrez.efetch("pubmed", id=number, retmode="xml")
    records = Entrez.parse(handle)
    a =  list(records)
    handle.close()
    return a


def summarise_record(records):
    "Short-cut function prints the title of any papers returned \
    by get_record()"
    for record in records:
        print(record['MedlineCitation']['Article']['ArticleTitle'])

        
# def expand(records):
#     if hasattr(records, '__iter__'):
#         if isinstance(records, list):
#             for item in records:
#                 expand(item)
#         else:
#             for item in records:
#                 expand(records[item])
#     else:
#         print records

def list_expand(records, st, dct):
    "list version of expand function - recursively calls itself and \
    dict version to fully investigate tree of list and dict items"
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
    """dict version of expand function - recursively calls itself and
    list version to fully investigate tree of list and dict items.
    """
    st = st + key + "\t"
    tdct = check(key, records, dct)
    #print dct
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
    """called by dict and list expand functions to compare."""
    if key in rdict.keys():
        rdict[key] = entry
    return rdict


def from_entrez(reference_dictionary, entrez_id):
    print "starting"
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
            "Title": "", # Journal title
            "Issue": "", # number
            "MedlinePgn": "",
}

# check("AuthorList", "hi there!")


#print "starting"
# with open('temp.txt', 'w') as f:  #
    # records = get_record("18716004")#,19304878,14630660")
    # st = ""
    # list_expand(records, st, f)
#from_entrez(ref_dict, "18716004")
#print "done"

#dict_convert(ref_dict, bib_dict)

#print bibtex_string("simon", bib_dict)
