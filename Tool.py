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
    #print key
    if hasattr(records, '__iter__'):
        if isinstance(records, list):
            for item in records:
                list_expand(item, st, handle)
        else:
            for item in records:
                dict_expand(item, records[item], st, handle)
    else:
        handle.write(st + "\t" + records +"\n")

def bibtex_string(bibtype, id, bibdict):
    output = bibtype + "{" + str(id) + ",\n"
    for key in bibdict:
        output += "\t" + key + " = {" + bibdict[key] + "},\n"
    return output + "}"


bibdict = {
    "author" : "",
    "address" : "",
    "annote" : "",
    "booktitle" : "",
    "chapter" : "",
    "crossref" : "",
    "edition" : "",
    "editor" : "",
    "howpublished" : "",
    "institution" : "",
    "journal" : "",
    "key" : "",
    "month" : "",
    "note" : "",
    "number" : "",
    "organization" : "",
    "pages" : "",
    "publisher" : "",
    "school" : "",
    "series" : "",
    "title" : "",
    "type" : "",
    "volume":"",
    "year":"",
}

print "starting"
with open('temp.txt', 'w') as f:
    records = get_record("18716004")#"19304878,14630660")
    st = ""
    list_expand(records,st,f)

print "done"
#summarise_record(records)

#full_info(records)

# count = 0
# for record in records:
#     print(record['MedlineCitation']['Article']['ArticleTitle'])
#     a = record
#     count += 1
