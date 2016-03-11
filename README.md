# Pubmed to Bibtex
A Python 2.7 tool that searches for a pubmed record and returns a bibtex citation.

Created as part of my MSc in Bioinformatics at Newcastle University

## Version History

**Version 0.0.1** (11/3/16) - current version

## Dependencies

This was built using Python 2.7 and [Biopython](http://biopython.org/wiki/Main_Page) 1.63 - use with other versions at your own risk

## Installation

Nothing too complex - simply unpack all three Python files into a single directory.

You will need to open pub_to_bibtex.py in an editor of your choice to set your email address.

The variable email is the first to be defined, and by default will be set to the string "email@email". If you do not change this the module will raise an exception.

**Some versions of this code exist with my email address "simon.c.baron@gmail.com" in this field - please do be sure to replace that with your own!**

Please abide by the [Entrez usage guildlines](http://www.ncbi.nlm.nih.gov/books/NBK25497/) when using this tool.

## Usage options

call the function using  Python from a terminal, eg:

**$**python pub_to_bibtex.py [arguments]

[arguments] to be replaced with options from among the following

Flag | Name | Description
------------ | ------------- | -------------
-h | help | Prints command line help showing options
[None] | id/search string | Any arguments not following a flag will be interpreted as pubmed IDs or an entrez search string
-f | file output | By default we display plain text in the command line, which can be copied into a .bib file, select this argument if you want to write to a text file.
--s | search | Enables searching functionality, replacing pubmed ids with a seach string
-n | number of entries | If searching is enabled, allows you to set the number of records to be returned from pubmed. Uses Entrez default of 20 if not specified.
-log | warnings | Create a log file for warnings generated
Any argument that comes directly after a -flag will be considered to be the argument of that flag, eg -f filename or -n number
Double flags (--s is the only one currently implemented) are options that do not require an arguement.

### Return citations by pubmed ID

eg:

    $python pub_to_bibtex.py 26944449

 return bibtex format for pubmed article with id 26944449

    $python pub_to_bibtex.py 26944449 18206625

sequentially return multiple article, ids 26944449 and 18206625

    python pub_to_bibtex.py 26944449 18206625 -f filename.bib

create (or overwrite) file "filename.bib" containing bibtex format citations for articles 26944449 and 18206625

### Search functionality

eg:

    $python pub_to_bibtex.py --s heatmaps

return 20 (default) citations about heatmaps

    $python pub_to_bibtex.py --s -n 3 heatmaps are pretty

return citations for the first 3 entrez search results for the phrase "heatmaps are pretty".

    $python pub_to_bibtex.py --s -n 3 heatmaps -f mybib.bib

create or overwrite a file (mybib.bib) containing the first 3 entrez search results for the phrase "heatmaps"

Search functionality is provided by Entrez esearch, details can be found on their [help page](http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch).

## Warnings

Warnings of type UserWarning are generated if a bibtex record does not meet expectations, eg an article which doesn't contain any title. These are printed to the shell. If you are having trouble reading them them I would recommend writing your bibtex to a file.

UserWarnings with message "failed to convert PMID:" indicate that one of the records has been unable to format correctly and will not be included.

This behaviour is designed to prevent a full crash if a record is not formatted in the expected manner.

-log filename can be used to create a file to store warning information.

## Contact

For help or bug reports or development suggestions, contact me at simon.c.baron@gmail.com.

## Future changes planned

- [x] Better error handling. As pubmed XML records are provided by authors who do not always follow the guidelines set down (and neither does the NCBI) error and warning handling are vital. Currently if you return 20 results and even one of them throws an error you get nothing, obviously this could be avoided.
- [x] Logging of warnings and errors generated
- [ ] Parsing an NCBI .txt download for collections or bibliograhies, allowing you to convert an existing PubMed collection to a bibtex bibliography.
- [ ] More complete text mining of PubMed XML - extend the information I can recover to more of bibtex fields.

Feel free to send me a branch request if you would like to help with any part of this.
