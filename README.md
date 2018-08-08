Bottom-up and top-down control of dispersal across major organismal groups
=====================================================================

Authors
---------------------------------------------------------------------
Emanuel A. Fronhofer, Delphine Legrand, Florian Altermatt, Armelle Ansart, Simon Blanchet, Dries Bonte, Alexis Chaine, Maxime Dahirel, Frederik De Laender, Jonathan De Raedt, Lucie di Gesu, Staffan Jacob, Oliver Kaltz, Estelle Laurent, Chelsea J. Little, Luc Madec, Florent Manzi, Stefano Masier, Felix Pellerin, Frank Pennekamp, Nicolas Schtickzelle, Lieven Therry, Alexandre Vong, Laurane Winandy and Julien Cote

Summary
---------------------------------------------------------------------
Data and code used within this [dispNet](https://dispnet.github.io/) project published under the title "Bottom-up and top-down control of dispersal across major organismal groups" in Nature Ecology and Evolution. A preprint version can be found on [bioRxiv](https://dx.doi.org/10.1101/213256).

Meta information
---------------------------------------------------------------------
The folder "data" contains the main data set as a text file. The "species" column specifies the species names, while RA and PRED refer to the treatment levels of resource availability and predation risk. Replicates and blocks are specified. "no_residents" and "no_dispersers" respectively refer to the measured number of residents (individuals remaining in the origin patch) and dispersers (individuals that have dispersed to the target patch). The data in the columns "relevant_taxon" and "lab" were used as random effects to correct for phylogeny and experimenters.

The folder "code" contains an R-script with a simulation of the model described in detail in the paper. This script reproduces the results shown in Fig. 2.