import os
from matplotlib import pyplot as plt
from scipy import signal
from scipy.stats import shapiro,bartlett,ttest_ind,wilcoxon,mannwhitneyu,kruskal,f_oneway
from Lib import getData,parse,DN_content,DN_skew,PullCodon,plot_stats,ParseAnnot,getData_aleatoire,analyze,DN_distribution,getLower,GC_skew_analysis,getRandom,Genome_signature_analysis
from Lib import GC_skew_statistical_analysis
import pandas as pd
import numpy as np

Species = ["Bacillius_anthracis","Oscilatoria_acuminata","Arabidopsis_thaliana"]
DNA = ["chr1","chr2","chr3","chr4","chr5","pOX1_bacilus_anthracis","pOSCIL1","pOSCIL2","Chlr","mito"]
labels = ["B.anthracis chromosome 1","B.anthracis pOX1","O.acuminata chromosome 1","O.acuminata pOSCIL1","O.acuminata pOSCIL2","A.thaliana chromosome 1","A.thaliana chromosome 2","A.thaliana chromosome 3","A.thaliana chromosome 4","A.thaliana chromosome 5","A.thaliana cpDNA","A.thaliana mtDNA"]
genomes_paths = [os.path.join(Species[0],DNA[0]+".txt"),os.path.join(Species[0],DNA[5]+".txt"),os.path.join(Species[1],DNA[0]+".txt"),os.path.join(Species[1],DNA[6]+".txt"),os.path.join(Species[1],DNA[7]+".txt"),os.path.join(Species[2],DNA[0]+".txt"),os.path.join(Species[2],DNA[1]+".txt"),os.path.join(Species[2],DNA[2]+".txt"),os.path.join(Species[2],DNA[3]+".txt"),os.path.join(Species[2],DNA[4]+".txt"),os.path.join(Species[2],DNA[8]+".txt"),os.path.join(Species[2],DNA[9]+".txt")]
random_seqs = ["10_BA_chr1","20_BA_chr1","40_BA_chr1","80_BA_chr1","10_BA_pOX1","20_BA_pOX1","40_BA_pOX1","80_BA_pOX1","10_OA_chr1","20_OA_chr1","40_OA_chr1","80_OA_chr1","10_OA_pOSCIL1","20_OA_pOSCIL1","40_OA_pOSCIL1","80_OA_pOSCIL1","10_OA_pOSCIL2","20_OA_pOSCIL2","40_OA_pOSCIL2","80_OA_pOSCIL2","10_BA_chr1","20_BA_chr1","40_BA_chr1","80_BA_chr1","10_BA_chr1","20_BA_chr1","40_BA_chr1","80_BA_chr1","10_BA_chr1","20_BA_chr1","40_BA_chr1","80_BA_chr1","10_BA_chr1","20_BA_chr1","40_BA_chr1","80_BA_chr1","10_BA_chr1","20_BA_chr1","40_BA_chr1","80_BA_chr1","10_AT_cpDNA","20_AT_cpDNA","40_AT_cpDNA","80_AT_cpDNA","10_AT_mtDNA","20_AT_mtDNA","40_AT_mtDNA","80_AT_mtDNA"]
random_seqs_paths = [os.path.join("Random_seqs",i+".txt") for i in random_seqs]

if __name__ == '__main__':

    ###############################
    ###### GC-skew analysis #######
    ###############################

    # Set window size to 0.2%:
    WS = 0.2
    # Analyze
    stds,means,UC_genomes,spectrums = GC_skew_analysis(genomes_paths,labels,WS,Species) #Pipelining Function (UC = Uppercased)
    GC_skew_statistical_analysis(spectrums, labels) #Proceed to statistical analysis with plot displaying


    #######################################
    ###### Genome signature analysis ######
    #######################################

    # Set the number of random sets used for each sequence
    n_random_sets = 4
    # Analyze
    Genome_signature_analysis(labels,random_seqs_paths,UC_genomes,n_random_sets)


































