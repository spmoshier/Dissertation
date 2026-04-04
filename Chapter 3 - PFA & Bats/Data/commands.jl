#R_HOME: /usr/local/R/gnu/14.1/4.4.0/lib64/R
#import Pkg; Pkg.add("RCall")
#Pkg.build("RCall")
#then using RCall loaded in no issue
#loaded in all R libraries within the online terminal

#module load Julia
#julia

using PhylogeneticFactorAnalysis

#check_beast()
#check_r()

pfa("simrun.xml")