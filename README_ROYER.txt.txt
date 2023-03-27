Below is the original header from the R code version of GEOCARBSULF.  This header file includes
several citations that went into the creation of GEOCARBSULF that will be useful in understanding
the equations used throughout the simulation.

Eric Mittelstaedt
March, 2023

# This r code is written by Dana Royer for the publication Royer et al (2014). The core GEOCARBSULF code (the second section of the code below) is translated mostly from the BASIC scripts written by Robert Berner in 2010 for the GEOCARBSULFvolc model (acquired by Dana Royer in October 2013), with updates from the FORTRAN scripts written by Jeffrey Park for Park & Royer (2011; see http://jparkcodes.blogspot.com/2011_06_01_archive.html); the most important component taken from the FORTRAN scripts is the convergence function for estimating CO2. 

#For r novices, running the model is straightforward. First, download the r program if you don't already have it (http://www.r-project.org/). Next, set the working directory (e.g., desktop) (File->Change dir...); this is where all files will be read from and written to. Next, put the two input csv files in your working directory. Finally, copy the entire content of the r-code file (i.e., the file you are reading right now) and paste it into the r console. The summary csv file will be outputted to your working directory automatically. A full, standard run (10000 resamples) takes ~5 minutes to run on a PC with a 3.1 GHz processor and 8 GB of RAM.

#There are two input files: the default values reproduce the light gray uncertainty envelopes in Royer et al (2014); to reproduce the dark gray envelopes, reduce all two sigma errors by 75%. "GEOCARB_input_summaries.csv" provides a summary of all input parameters (see also Appendices 1-3 in Royer et al, 2014), including all necessary information for the time-invariant constants. "GEOCARB_input_arrays.csv" provides the time-step mean and two sigma values for the subset of parameters that are time arrays. Please note: if you open GEOCARB_input_summaries.csv in Excel, the "-inf" cells will probably be read in as an error; to correct, type <<'-inf>>.
#In the "input summaries" file:
#  parameters ending in "_Godderis" are only activated if "Godderis" is set to TRUE (third line in the code below); in this case, the corresponding non-Godd?ris parameter choices are ignored.
#  the "type" column differentiates time-invariant parameters (constants) from time-dependent parameters (arrays).
#  the "resample" column describes whether a particular parameter will undergo resampling for error analysis. If you wish to forgo all resampling, set resampleN to 1 (first line in the code below); this will override whatever is written in this particular column.
#  the "distribution" column describes whether the resampled distributions will follow a normal or lognormal distribution.
#  the "lower_limit" and "upper_limit" columns describe whether resamples should be clipped at particular threshold values; in these cases, a value just inside the limit (+/- 0.0001) is assigned; in the default model run, this rarely happens.

#Useful references:
#Berner RA. 2004. The Phanerozoic Carbon Cycle. Oxford University Press.
#Berner RA. 2006a. GEOCARBSULF: a combined model for Phanerozoic atmospheric O2 and CO2. Geochimica et Cosmochimica Acta 70: 5653-5664.
#Berner RA. 2006b. Inclusion of the weathering of volcanic rocks in the GEOCARBSULF model. American Journal of Science 306: 295-302.
#Berner RA. 2008. Addendum to "Inclusion of the weathering of volcanic rocks in the GEOCARBSULF model". American Journal of Science 308: 100-103.
#Godd?ris Y, Donnadieu Y, Lefebvre V, Le Hir G & Nardin E. 2012. Tectonic control of continental weathering, atmospheric CO2, and climate over Phanerozoic times. Comptes Rendus Geoscience 344: 652-662.
#Park J & Royer DL. 2011. Geologic constraints on the glacial amplification of Phanerozoic climate sensitivity. American Journal of Science 311: 1-26.
#Royer DL, Donnadieu Y, Park J, Kowalczyk J, Godd?ris Y. 2014. Error analysis of CO2 and O2 estimates from the long-term geochemical model GEOCARBSULF. American Journal of Science, 314: 1259-1283..