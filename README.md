# beltman-DRAGbarcode-cleanup

Contains:

-Scripts to remove spurious barcodes in barcoding deepseq data generated from the genetic 'barcode mouse' construct as generated in the lab of Ton Schumacher (NKI, NL)

-Example input and output files


Usage:

-Within the script, make a correct reference to the files you aim to clean up.

-Place the files to be cleaned in the 'rawdata' directory, and adjust the file 'samplenames.txt' (also in that directory). For details on what is expected see the remarks in the header of 'tagcleanup.R' and the example files.

-Run the script 'tagcleanup.R'.

-The output files are placed in the directory 'output' (for details see the header of 'tagcleanup.R' and the example output that results from running the script on the example input).

This approach to clean DRAG barcoding data is used in https://www.biorxiv.org/content/10.1101/2022.12.06.519273v1
