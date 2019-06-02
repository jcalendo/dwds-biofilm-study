# dwds-biofilm-study

This is the project directory for all of the work I have completed so far.

## /bin

This directory contains all scripts and notebooks needed to reproduce the analysis in its entirety. It contains three subdirectories:

 - html_notebooks/ is a directory containing html versions of the jupyter notebooks that can be opened in a web browser
 - notebooks/ contain jupyter notebooks with functional code. The code in these notebooks is what performed the bulk of the analysis. They also incldue some explanations of each analysis step. If you wish to run these scripts yourself, be mindful of your current working directory.
- scripts/ contains R scripts used for analysis and plotting

## data/

This directory contains the raw reads. This is the basis for all later steps and the most important directory. It is split into two subdirectories:

 - batch_1_reads/ the first sequencing run
 - batch_2_reads/ the second sequencing run

## doc/

The /doc fold contains all text documents relating to the project. It has four subdirectories:

 - final-paper/ contains all a pdf of my final paper as well as a folder containing all LaTex code and images used to produce the paper.
 - final-presentation/ contains a pdf version of my final capstone presentation
 - notes/ contains miscellaneous meeting notes
 - paperwork/ contains all relevant capstone paperwork

## results/ 

The results directory contains all results produced at every stage of the experiment. **The latest results are located in 2019-04-01**. Use these results! Each of the results folders contains many subdirectories. The meaning of the files and directories in these folders is explained in the notebooks files that are used to create the results. 

## src/

Currently empty. The plan was to create a master script that could re-run the entire analysis - that script would have been located here. 
