"""Command line module for converting pubmed to bibtex.

Requires Fetch_Record and Format_Record in the same directory.
"""

import Fetch_Record
import Format_Record

import argparse

parser = argparse.ArgumentParser(description='Get a pubmed file and print a bibtex format citation.')

parser.add_argument('pubmed_id', metavar='id/search string',
                    type=str, nargs='+', help='pubmed id to be retrieved')


parser.add_argument('--s', dest='search',
                    action='store_const', const=True, default=False,
                    help='search Entrez for a string '
                    '(Default: fetch specific ID)')

parser.add_argument('-n', dest='n', , nargs=1,
                    help="number of records if searching")

args = parser.parse_args()



if args.search:
    entrez_input = " ".join(args.pubmed_id)
else:
    entrez_input = ",".join(args.pubmed_id)

print entrez_input

print args.n

#Format_Record.format_convert(entrez_input, search=args.search, retmax=args.n)


# to file?
