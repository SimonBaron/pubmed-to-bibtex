'''testings some thgins to see if I can make a tool'''

from Bio import Entrez

Entrez.email = "simon.c.baron@gmail.com"

handle = Entrez.efetch("pubmed", id="19304878,14630660", retmode="xml")
records = Entrez.parse(handle)


count = 0
for record in records:
    print(record['MedlineCitation']['Article']['ArticleTitle'])
    a = record
    count += 1

    
print "hello"

handle.close()
