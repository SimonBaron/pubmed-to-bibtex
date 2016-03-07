import Fetch_Record as fetch

def bibtex_string(id, bib_dict):
    output = "@" + bib_dict["type"] + "{" + str(id) + ",\n"
    for key in bib_dict:
        output += "\t" + key + " = {" + bib_dict[key] + "},\n"
    return output + "}"

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

citations = fetch.from_entrez(fetch.ref_dict, "18716004")

dict_convert(citations[0], bib_dict)

print bibtex_string("yay", bib_dict)
