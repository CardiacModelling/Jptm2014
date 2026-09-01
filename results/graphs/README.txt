This is a set of instructions for generating the results shown in the paper.

process_results.m
 - loads up the raw C++ output,
 - plots the action potential duration against concentration curves
 - calculates the contingency table entries.
 - stores the contingency table 'tally' in a file called results.mat

plot_results.py 
 - matplotlib was used in the end to make the pretty figures in the paper
   (since Matlab proved a bit useless at it). Thanks to Tom Dunton for this!

results_to_latex_tables.m
 - loads up results.mat
 - generates latex format output for pasting into the paper.

input_data_to_latex_table.m
 - loads up the input data from the folder in the Jptm2014 project.
 - generates latex format output for pasting into the paper.
