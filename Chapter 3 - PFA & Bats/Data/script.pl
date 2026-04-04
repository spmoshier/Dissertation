#!/usr/bin/perl

 

use POSIX;

 

print "\n\n\n Simulation starting ..... \n\n\n";

 

# In the .xml in this directly, we need the trait data in file simdat1.csv

# and the tree data in simtree1.csv

# trait data: simdat1.csv

# tree files: simtree1.txt

 

# Also need a file called commands.jl that contains the julia commands to run pfa:

# using PhylogeneticFactorAnalysis

# pfa("Julia_runs/julia_sim/simrun.xml")

 

 

use constant NLOOP1 => 2402;  

$count = 1;

 

while ($count < NLOOP1)

        {

 

                system("cp simdat$count.csv simdat.csv");

                system("cp simtree$count.txt simtree.txt");

                system("julia < commands.jl > julia_output.txt");

		system("cp julia_output.txt pfa_results/julia_output$count.txt");

                system("cp pfa_output/pfa_output_loadingsStatistics.csv pfa_results/loadingsStats$count.csv");

                system("cp pfa_output/pfa_output_factorMeans.csv pfa_results/factorMeans$count.csv");

		system("cp pfa_output/pfa_output_loadings.pdf pfa_results/loadings$count.pdf");

 

                $count = $count+1;

 

}

 

print "Done!\n\n";

 

#######################################################

# Examples below for calling R and julia:

#system("R CMD BATCH testscript.R");

#system("julia < testcommands.jl > juliaoutput.txt");

######################################################