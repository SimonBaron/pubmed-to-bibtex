"""Formatting module for pubmed>bibtex converter.

Contains the functions necessary to convert a dictionary provided by the
Fetch_Record module and format them in a bibtex friendly format.

Globally defined dictionaries bib_dict and pubtype_dict act as templates.

bib_dict is a dictionary which when populated can be printed as a bibtex
string,

pubtype_dict contains the conversions of pubmed record types into bibtex
types.

Includes basic warning functionality.
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
    """Convert data from pubmed-like dictionary into bibtex-like dict.

    calls auxilliary functions for more complex conversions: type and name.

    pubmed_info is a dict of fetch.ref_dict format, containing the parsed
    information from a pubmed record.

    bib_dict is a COPY of globally defined bib_dict which will be populated
    and returned.

    type_ref_dict is a dictionary of type conversions such as pubtype_dict.

    known issue - if author has not included "ForeName" will crash.
    """
    new_dict = bib_dict
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
    new_dict["author"] = assign_author(pubmed_info["AuthorList"],
                                       pubmed_info["PMID"])
    new_dict["type"] = assign_type(type_ref_dict,
                                   pubmed_info["PublicationTypeList"])
    return new_dict


def assign_author(pubauthor_list, pmid):
    """Format a string of author names.

    Pubmed type is to have a list of dictionaries, one
    for each author. Bibtex format takes a single string."""
    author_string = ""
    for author in pubauthor_list:
        try:
            FN = author["ForeName"]
        except KeyError:
            error = "{} has an author with no forename".format(pmid)
            warnings.warn(error)
            FN = ""
        try:
            LN = author["LastName"]
        except KeyError:
            error = "{} has an author with no lastname".format(pmid)
            warnings.warn(error)
            LN = ""
        author_string += u"{} {} and ".format(FN, LN)
    return author_string[:-5]   # slice off the final " and "

def assign_type(pubtype_dict, pubtype_list):
    """Assign a bibtex type to a list of  pubmed types."""
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


def format_convert(entrez_key, email, search=False, retmax=3):
    """Converts parametric input to functions.

    passes inputs for entrez term (id or search strings),
    email (string format email of querier),
    search boolean selection and retmax to correct functions to output
    a list of formatted bibtex strings.

    try and except used to downgrade errors for a single record into
    warnings for the entire process.
    """
    citations = fetch.from_entrez(fetch.ref_dict,
                                  entrez_key, email, retmax, search)
    output = []
    for record in citations:
        try:
            c_record = dict_convert(record, bib_dict.copy(), pubtype_dict)
        except:
            message = "failed to convert PMID: {}".format(record["PMID"])
            message += "\ncontinuing with other records."
            warnings.warn(message)
        output.append(bibtex_string(c_record))
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
