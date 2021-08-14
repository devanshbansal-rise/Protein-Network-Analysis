# Protein-Network-Analysis
Tool for Computational Protein Network Analysis

This file contains code that can be run where R code is run to create an interactive app that computes various forms of computational protein network analysis, some of which include computing comparison and overlap, computing differential connectivity of proteins, identifying key proteins, enriching biological pathways of key proteins, and visualizing networks. 

The interactivity of the app allows it to accept any user data and generate networks and analysis. The data to be input must be in datatable/dataframe format, each datatable representing a network and formatted the same way with the same protein IDs. That format must contain two columns with the interacting proteins, each row of those two columns representing an interaction in the network. Inputting just two networks may result in some errors in the comparison/overlap section of analysis, and inputting more than 6-10 networks will make the app very slow; the recommended number of networks to input is around 3-6. This analysis only works for unweighted and non-directional networks, and will not consider weights or directionality when the data is input.

ID mapping is conducted in this tool just as it is on the UniProt ID Mapping website, except the tool automatically maps to the UniProtKB-ID gene names for the main purpose of enrichment through EnrichR. The mapping can be done within the program itself using the mapping table downloaded from UniProt's website, which you must download yourself from https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ (the last file titled "idmapping_selected_tab) and rename to "mapping.tab" for the app to read. ID mapping is only necessary for biological pathway enrichment AND only if your IDs are not already gene names, so this step isn't necessary if you don't need it.

This repository also contains some of the human tissue interactomes found in the TissueNet database that can be uesd to test the tool (they use columns 1 and 2 of the datatables and map from Ensembl protein IDs).

Descriptions about each part of the analysis and any remaining information necessary can be found within the app itself.
