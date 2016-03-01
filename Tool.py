'''testings some things to see if I can make a tool'''

from Bio import Entrez

Entrez.email = "simon.c.baron@gmail.com"

def get_record(number):

    handle = Entrez.efetch("pubmed", id=number, retmode="xml")
    records = Entrez.parse(handle)
    a = list(records)
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


print "starting"

records = get_record("14716003")#"19304878,14630660")

summarise_record(records)

full_info(records)

# count = 0
# for record in records:
#     print(record['MedlineCitation']['Article']['ArticleTitle'])
#     a = record
#     count += 1




