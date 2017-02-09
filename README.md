# EBI-Metagenomics-collector
Author: Mike Aaldering
Company: NIOO Knaw
Email: M.Aaldering2@nioo.knaw.nl
Private Email: mike.aaldering098@gmail.com
Download all metadata and OTU tables from metagenome samples.


To run EBI.py , the following python libraries are needed: 

-xml 
-biopython
-pandas 
-unicodedata 
-urllib2 


The following .txt files are needed:

-"InterestingSpecies.txt" Type here the species that you need, all the other species will be filtered out.

-"Metadata_order.txt" Make a list with metadata names and the columns in the .csv output file will have the same order.




Search metadata and species from EBI metagenomics and store it to a .csv file.
1       Download metadata on ENA.
    1.1 Download run, sample and project metadata on EBI metagenomics.
2       Search article id's for every sample on pubmed central and add article id to metadata.
3       Download UTO table for every sample.
4       Save  the metadata and OTU table to a seperate CSV file.
5       Normalize metadata.
6       Join OTU table and metadata and save to .csv file.
