"""Command line module for converting pubmed to bibtex.

Requires Fetch_Record and Format_Record in the same directory.

Basic help available, use:
$python pub_to_bibtex.py -h

For more usage and documentation please see the readme file available from
https://github.com/SimonBaron/pubmed-to-bibtex

Please remember to replace the email variable with your correct email address
This is passed to Entrez along with queries.
"""

import Fetch_Record

import Format_Record

import argparse

import logging

class EmailError(Exception):
    """Raise an error if email is not set correctly."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

email = "email@email"

if email == "email@email":
    raise EmailError("You should always tell Entrez who you are - "
                     "please edit pub_to_bibtex.py and "
                     "replace the string with your own")

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

parser.add_argument('-f', metavar='file_name', type=str, nargs='?',
                    help="if present will write output to file")

parser.add_argument('-log', metavar='log file', type=str, nargs='?',
                    help="if present will write logs to file")

args = parser.parse_args()

if args.log is not None:
    logging.basicConfig(filename=args.log)
    logging.captureWarnings(True)
    print "writing logs to {}".format(args.log)

if args.search:
    entrez_input = " ".join(args.pubmed_id)
else:
    entrez_input = ",".join(args.pubmed_id)


bibtex = Format_Record.format_convert(entrez_input, email,
                                      search=args.search, retmax=args.n)

if args.f is not None:          # if -f filename selected
    print "writing output to file: {}".format(args.f)
    with open(args.f, 'w') as myfile:
        for record in bibtex:
            myfile.write(record.encode('utf8'))
else:
    for record in bibtex:
        print record
