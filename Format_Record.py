"""Formatting module for pubmed>bibtex converter.

Calls the Fetch_Record module to get pubmed records via entrez,
and formats them in a bibtex friendly format.

Currently only supports article type documents.
"""

import Fetch_Record as fetch


def bibtex_string(id, bib_dict):
    """Reformat a dictionary into a string of bibtex format."""
    output = "@" + bib_dict["type"] + "{" + str(id) + ",\n"
    for key in bib_dict:
        output += "\t" + key + " = {" + bib_dict[key] + "},\n"
    return output + "}"


def dict_convert(pubmed_info, bib_dict, type_ref_dict):
    """Convert data from pubmed-like dictionary into bibtex-like dict."""
    new_dict = bib_dict
    author_string = ""
    print pubmed_info["AuthorList"]
    for i in range(len(pubmed_info["AuthorList"])):
        author_string += "{} {} and ".format(
            pubmed_info["AuthorList"][i]["ForeName"],
            pubmed_info["AuthorList"][i]["LastName"])
    new_dict["author"] = author_string[:-4]
    new_dict["month"] = pubmed_info["PubDate"]["Month"]
    new_dict["year"] = pubmed_info["PubDate"]["Year"]
    new_dict["volume"] = pubmed_info["Volume"]
    new_dict["title"] = pubmed_info["ArticleTitle"]
    new_dict["journal"] = pubmed_info["Title"]
    new_dict["number"] = pubmed_info["Issue"]
    new_dict["type"] = assign_type(type_ref_dict,
                                   pubmed_info["PublicationTypeList"])
    return bib_dict


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
    "organization": "",
    "pages": "",
    "publisher": "",
    "school": "",
    "series": "",
    "title": "",
    "type": "",
    "volume": "",
    "year": "",
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

citations = fetch.from_entrez(fetch.ref_dict, "heatmap", search=True, retmax=3)

for citation in citations:
    populated_dict = dict_convert(citation, bib_dict.copy(), pubtype_dict)
    print bibtex_string("yay", populated_dict)
