This git repository contains some sample data files and a Makefile for running the python (v2.7) script.

Two examples are included and can be generated as follows: 
- **make plot** should generate some sample CircleMap image files into a folder called Circle output.
- **make -f circleMapWithColorFile.mak plot** should generate CircleMap image files into the Circle output folder using custom colors (specified in thefile Data/kegg-slea.colorcoding.tsv)

To plot your own CircleMaps, you can put your own data files in the data directory and then make the corresponding changes to the makefile.

To clean the directory (remove results and temporary files generated upon execution) type the following command: **make clean**

CITATION
Wong CK, Vaske CJ, Ng S, Sanborn J, Benz S. Haussler D, Stuart J. The UCSC Interaction Browser: Multi-dimensional data views in pathway context. Nucleic Acids Research. 2013 Jul 1;41(Web Server issue):W218-24. doi: 10.1093/nar/gkt473. PMID:23748957.
