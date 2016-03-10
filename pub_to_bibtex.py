"""Command line module for converting pubmed to bibtex.

Requires Fetch_Record and Format_Record in the same directory.

parser design:
1. get record given id
 - should support multiple ids
 - should support writing to file

2. search given a string
 - needs to return n records

eg:

1. python pub_to_bibtex.py 239487 987348 987343
   python pub_to_bibtex.py -f mybib.bib 239487 987348 987343

2. python pub_to_bibtex.py -s "the title of some paper" (-f)

test input = 26944449 18206625

"""

import Fetch_Record

import Format_Record

import argparse

parser = argparse.ArgumentParser(
    description='Get a pubmed file and print a bibtex format citation.')

parser.add_argument('pubmed_id', metavar='id/search string',
                    type=str, nargs='+', help='pubmed id to be retrieved')


parser.add_argument('--s', dest='search',
                    action='store_const', const=True, default=False,
                    help='search Entrez for a string '
                    '(Default: fetch specific ID)')

parser.add_argument('-n', metavar='n', type=int, nargs='?',
                    help="number of records(only used if searching)")

parser.add_argument('--f', metavar='file_name', type=str, nargs='?',
                    help="if present will write to file")

args = parser.parse_args()

if args.search:
    entrez_input = " ".join(args.pubmed_id)
else:
    entrez_input = ",".join(args.pubmed_id)

# print entrez_input

# print args.n

# print args.f is not None

bibtex = Format_Record.format_convert(entrez_input,
                                      search=args.search, retmax=args.n)

if args.f is not None:
    print "writing to file: {}".format(args.f)
    with open(args.f, 'w') as myfile:
        for record in bibtex:
            myfile.write(record.encode('utf8'))
else:
    for record in bibtex:
        print record
