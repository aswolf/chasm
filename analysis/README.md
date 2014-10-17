#analysis directory
This directory contains subdirs each representing a certain type of analysis.
The analysis subdirs (analysis-1, analysis-2, etc...) contain scripts that are run to produce figures, tables, etc.
Each subdir has its own figs and tbls directories to organize output.
Most of these steps are exploratory data analysis, and will be contained in a collection of scripts, which are all logically named.
The best-case scenario is that the analysis is contained in an ipython notebook (.ipynb file).
If so, the ipynb file must be cleaned of output before storing in a git repo.

#version control
For space reasons, the exploratory figures are not kept under version control.
Tables are small, so they will be.
If version control is desired for figures, they should be copied to the writeup directory, where the manuscript is stored.
The manuscript figures are version controlled, since they are quite important, and the manuscript should be buildable from the repo exactly.
