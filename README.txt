Chaste "Jptm2014" bolt-on project to accompany the paper:
"Prediction of Thorough QT trial results using action potential simulations based on ion channel screens"

This Chaste project is written to work with version 3.2 of Chaste only.

To use this project you now need to checkout the ApPredict project too.

Make sure both projects are either present in Chaste/projects or have symbolic links there.

INPUT DATA
The three high throughput screening datasets that were used in the paper
AZ_HTS_data.txt - 'Q' dataset,
GSK_HTS_data.txt - 'B&Q2' dataset,
AZ_HTS_data_Gintant_hERG.txt - 'M&Q' dataset.

SIMULATION CODE
This project contains just two files for running the simulations 
(the rest are pulled in from the public ApPredict project). Note that
the ApPredict project is kept up to date with the development version of Chaste
so you need to obtain the Chaste 3.2 version of ApPredict too).

test/TestMakeLookupTable.hpp - this pre-computes the APD for as many combinations of
the 5 channel blocks as you give it chance to (a handful of settings to change in the
first few lines - model, pacing rate, number of parameters to sweep).
To get to around a million entries, which gives good coverage in 5D takes around
a month on a 12 core machine.
Note that this test checkpoints, and so can be stopped and restarted (just make sure 
it isn't in the middle of writing a checkpoint when you stop it!).
This test can be ignored, and if so you will get predictions without credible regions
when running the main simulation test.

test/TestTqtCompounds - this test makes an executable that can be run with a variety
of options to recreate the results shown in the paper. 

RESULTS
We have included the raw results of the simulation study, so they can be re-analysed without 
re-running all the simulations if you wish. A couple of Matlab scripts do this, another readme
is in the results/graphs folder.

The TQT results are in the text file 
results/graphs/clinical_tqt_results 
these are taken from the Supplementary Material S2.


  
