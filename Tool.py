'''testings some things to see if I can make a tool'''

from Bio import Entrez

Entrez.email = "simon.c.baron@gmail.com"


def get_record(number):

    handle = Entrez.efetch("pubmed", id=number, retmode="xml")
    records = Entrez.parse(handle)
    a =  list(records)
    handle.close()
    return a


def summarise_record(records):
    for record in records:
        print(record['MedlineCitation']['Article']['ArticleTitle'])

        
def full_info(records):
    for record in records:
        for key in record["MedlineCitation"].keys():
            print key,
            print "\t"
            print record["MedlineCitation"][key]

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

def list_expand(records, st,handle):
    if hasattr(records, '__iter__'):
        if isinstance(records, list):
            for item in records:
                list_expand(item, st, handle)
        else:
            for item in records:
                dict_expand(item, records[item], st, handle)
    else:
        handle.write(st + records+"\n")

def dict_expand(key, records, st, handle):
    st = st + key + "\t"
    check(key, records)
    if hasattr(records, '__iter__'):
        if isinstance(records, list):
            for item in records:
                list_expand(item, st, handle)
        else:
            for item in records:
                dict_expand(item, records[item], st, handle)
    else:
        handle.write(st + "\t" + records + "\n")

        
def bibtex_string(id, bib_dict):
    output = "@" + bib_dict["type"] + "{" + str(id) + ",\n"
    for key in bib_dict:
        output += "\t" + key + " = {" + bib_dict[key] + "},\n"
    return output + "}"


def check(key, entry):
    global ref_dict
    if key in ref_dict:
        ref_dict[key] = entry

def dict_convert(ref_dict, bib_dict):
    author_string = ""
    for i in range(len(ref_dict["AuthorList"])):
        author_string += "{} {} and ".format(
            ref_dict["AuthorList"][i]["ForeName"],
            ref_dict["AuthorList"][i]["LastName"])
    bib_dict["author"] =  author_string[:-4]
    bib_dict["month"] = ref_dict["DateCompleted"]["Month"]
    bib_dict["year"] = ref_dict["DateCompleted"]["Year"]
    bib_dict["volume"] = ref_dict["Volume"]
    bib_dict["title"] = ref_dict["ArticleTitle"]
    bib_dict["journal"] = ref_dict["Title"]
    bib_dict["number"] = ref_dict["Issue"]
    if str(ref_dict["PublicationTypeList"][0]) == "Journal Article":
        bib_dict["type"] = "article"






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

bib_dict = {
    "type": "", ##
    "author": "", ###
    "address": "",
    "annote": "",
    "booktitle": "",
    "chapter": "",
    "crossref": "",
    "edition": "",
    "editor": "",
    "howpublished": "",
    "institution": "",
    "journal": "", ##
    "key": "",
    "month": "", ###
    "note": "",
    "number": "", ##
    "organization": "",
    "pages": "", ##
    "publisher": "",
    "school": "",
    "series": "",
    "title": "", ##
    "type": "",
    "volume": "", ###
    "year": "", ###
}

print "starting"
with open('temp.txt', 'w') as f:
    records = get_record("18716004")  # "19304878,14630660")
    st = ""
    list_expand(records, st, f)

print "done"

dict_convert(ref_dict, bib_dict)

print bibtex_string("simon", bib_dict)
