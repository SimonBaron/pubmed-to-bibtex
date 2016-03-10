"""Formatting module for pubmed>bibtex converter.

Contains the functions necessary to convert a dictionary provided by the
Fetch_Record module and format them in a bibtex friendly format.
"""

import Fetch_Record as fetch

import warnings


def bibtex_string(bib_dict, id="default"):
    """Reformat a dictionary into a string of bibtex format."""
    check_quality(bib_dict)
    if id == "default":
        id = bib_dict["id"]
    output = "@" + bib_dict["type"] + "{" + str(id) + ",\n"
    for key in bib_dict:
        if bib_dict[key] != "":
            output += "\t" + key + " = {" + bib_dict[key] + "},\n"
    return output + "}\n\n"


def dict_convert(pubmed_info, bib_dict, type_ref_dict):
    """Convert data from pubmed-like dictionary into bibtex-like dict."""
    new_dict = bib_dict
    author_string = ""
    for i in range(len(pubmed_info["AuthorList"])):
        author_string += u"{} {} and ".format(
            pubmed_info["AuthorList"][i]["ForeName"],
            pubmed_info["AuthorList"][i]["LastName"])
    new_dict["author"] = author_string[:-5]
    if "Month" in pubmed_info["PubDate"]:
        new_dict["month"] = pubmed_info["PubDate"]["Month"]
    if "Year" in pubmed_info["PubDate"]:
        new_dict["year"] = pubmed_info["PubDate"]["Year"]
    new_dict["volume"] = pubmed_info["Volume"]
    new_dict["title"] = pubmed_info["ArticleTitle"]
    new_dict["journal"] = pubmed_info["Title"]
    new_dict["id"] = pubmed_info["PMID"]
    new_dict["number"] = pubmed_info["Issue"]
    new_dict["abstract"] = pubmed_info["AbstractText"][0]
    new_dict["type"] = assign_type(type_ref_dict,
                                   pubmed_info["PublicationTypeList"])
    return new_dict


def assign_type(pubtype_dict, pubtype_list):
    """Assign a bibtex type to a pubmed typelist."""
    possibles = []
    for pubtype in pubtype_list:
        if pubtype.lower() in [word.lower() for word in pubtype_dict.keys()]:
            possibles.append(pubtype_dict[pubtype.lower()])
    if "article" in possibles:
        return "article"
    if "inproceedings" in possibles:
        return "inproceedings"
    if "proceedings" in possibles:
        return "proceedings"
    if "misc" in possibles:
        return "misc"
    else:
        return "nonsense"


def check_quality(mydict):
    """Check the contents of a bibtex-type dictionary before printing it.

    Bibtex format will not crash if used with insufficient data, but does
    expect certain fields to be populated depending on the type of record.

    Raises UserWarning giving the type of the record and the missing field.
    """
    contains = []
    for entry in mydict:
        if mydict[entry] != "":
            contains.append(entry)
    if mydict["type"] == "article":
        required = ["year", "title", "journal", "author", "volume"]
    if mydict["type"] == "inproceedings":
        required = ["year", "title", "booktitle", "author"]
    if mydict["type"] == "proceedings":
        required = ["year", "title"]
    if mydict["type"] == "misc":
        required = []
    if mydict["type"] == "nonsense":
        required = ["year", "title", "journal", "author", "volume"]
    for entry in required:
        if entry not in contains:
            message = ("In record with pubmed id {}: \n{}-type record "
                       "would expect to contain an instance of {}").format(
                           mydict["id"], mydict["type"], entry)
            warnings.warn(message)

def format_convert(entrez_key, search=False, retmax=3):
        citations = fetch.from_entrez(fetch.ref_dict,
                                      entrez_key, retmax, search)
        output = []
        for record in citations:
             output.append(bibtex_string(
                dict_convert(record, bib_dict.copy(), pubtype_dict)))
        return output

bib_dict = {
    "type": "",
    "author": "",
    "address": "",
    "annote": "",
    "booktitle": "",
    "chapter": "",
    "crossref": "",
    "edition": "",
    "editor": "",
    "howpublished": "",
    "institution": "",
    "journal": "",
    "key": "",
    "month": "",
    "note": "",
    "number": "",
    "organization": "",         # could be done - who hosted the conference?
    "pages": "",
    "publisher": "",
    "school": "",
    "series": "",
    "title": "",
    "type": "",
    "volume": "",
    "year": "",
    "id": "",
    "abstract": ""
    }

pubtype_dict = {
    "addresses": "proceedings",
    "bibliography": "misc",
    "case reports": "misc",
    "classical article": "article",
    "clinical conference": "inproceedings",
    "clinical trial": "misc",
    "congresses": "inproceedings",
    "consensus development conference": "inproceedings",
    "consensus development conference, NIH": "inproceedings",
    "corrected and republished article": "article",
    "editorial": "article",
    "festschrift": "article",
    "guideline": "misc",
    "interview": "misc",
    "journal article": "article",
    "lectures": "misc",
    "letter": "misc",
    "meta-analysis": "article",
    "news": "misc",
    "newspaper article": "article",
    "observational study": "misc",
    "patient education handout": "misc",
    "practice guideline": "misc",
    "published erratum": "article",
    "retraction of publication": "article",
    "review": "article",
    "video-audio media": "misc",
    "webcasts": "misc"
    }

if __name__ == "__main__":
    citations = fetch.from_entrez(fetch.ref_dict,
                                  "heatmap", search=True, retmax=3)
    for citation in citations:
        populated_dict = dict_convert(citation, bib_dict.copy(), pubtype_dict)
        print bibtex_string(populated_dict)
