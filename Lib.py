#Count(G + C)/Count(A + T + G + C) * 100 = GC content
#GC_skew = Count(G-C)/G+C

from scipy.stats import shapiro,bartlett
from scipy.stats import ttest_ind,wilcoxon,mannwhitneyu,kruskal,f_oneway
from scipy import stats
import os
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from scipy.fftpack import fft
from scipy import signal
import numpy as np


# Pull datas from text file
def getData(path):
    # Open the texte file
    genome_file = open(path, "r")
    # Get file_size
    genomeSize = os.stat(path).st_size
    # Get file lines
    lines = genome_file.readlines()
    lines_list = [line for line in lines]
    # Remove the first line
    lines_list.pop(0)
    return lines_list,genomeSize


# Define a parsing function
def parse(lines_list):
    bases = ["A","T","G","C"]
    genome = [i for j in lines_list for i in j if i in bases] #Pull characters from line_list on if there are in bases (ie : are nucleotides)

    return genome


# Define a DN_content getter
def DN_content(genome,label):
    bases = ["A","T","G","C"]
    counts = [genome.count(i) for i in bases] # Count each nucleotides in the genome
    print("##################################### : ",label)
    for i in range(len(counts)):
        print(bases[i]," Counts : ",counts[i], "   % : ",round(counts[i]/sum(counts)*100)) # Compute their distributions
    AT_content = (counts[0] + counts[1]) / sum(counts) * 100
    GC_content = (counts[2]+counts[3])/sum(counts)*100

    print("GC content : ", round(GC_content),"%")
    print("AT content : ", round(AT_content),"%")


# Define a DN_skew getter
def DN_skew(genome,fraction):
    #we define an inital window with a size of 0.2%
    window_size = round((len(genome)*fraction)/100)
    skews = []
    position = []
    window = genome[0: window_size]
    #Define GC-skew formula
    nG = window.count("G")
    nC = window.count("C")
    cumulative = []
    GC_skew = (nG - nC) / (nG + nC)
    cumulative.append(GC_skew)

    # Calculate GC-skews along the sequence using the defined window
    for i in range(1,len(genome)-window_size,window_size):
        window = genome[i:i+window_size]
        nG = window.count("G")
        nC = window.count("C")

        try: # try to calculate GC_skew
            GC_skew = (nG-nC)/(nG+nC)
            cumulative.append(GC_skew+cumulative[-1])
        except: # Pass if nG==nC (no skews)
            ZeroDivisionError
        skews.append(GC_skew)
        position.append(i)
    return skews,position,cumulative


# Function to plot GC_skews
def plot_skew(genome,specie,idx,WS):
    ###Generate the lagging strand

    dic = {"A": "T", "T": "A", "G": "C", "C": "G"}
    #Get complementary
    reverse = [dic.get(i) for i in genome if i in dic]
    #Get Reverse
    reverse.reverse()
    #Calculate GC_skew
    a2, b2,c2 = DN_skew(genome, WS) #brin sens
    a4, b4, c4 = DN_skew(reverse, WS)  # brin antisent

    ###Apply fast fourrier transform to obtain (not applied here)
    #a1t= fft(c1) #Fast fourrier transform of cumulative GC-skew
    #a2t = fft(c3) #Fast fourrier transform of cumulative GC-skew

    ###Display annotations : Pull gene position (not applied here)
    #gene_positions = [int(i) for i in gene_position if i[-1] != '.' and i[-1] != ',']
    #print("La position des genes est : ", gene_positions)

    ### Compute some stats
    std = np.std(a2) #std
    abs_a2 = [abs(i) for i in a2]#absolute stf
    mean = np.mean(abs_a2)#mean

    ###Create the plot
    plt.figure(1,figsize=(16,16))
    plt.suptitle(specie, fontsize="x-large")

    ##First panel : Leading strand GC_skew
    plt.subplot(221)
    plt.plot(b2,a2,'r-')
    plt.title("GC-skew")
    plt.ylabel("Leading strand : WS = 0.2 %")
    ##Display gene position (not applied here)
    # plt.scatter(gene_positions, [i * 0 for i in range(len(gene_positions))],c='red',s=0.25)

    ##Second panel : Leading strand cumulative GC-skew
    plt.subplot(222)
    plt.plot(range(len(c2)), c2)
    plt.title("Cumulative GC-skew")
    plt.ylabel(" WS = 0.2%")
    ##Third panel : Lagging strand GC_skew
    plt.subplot(223)
    plt.plot(b4,a4,'g-')
    plt.xlabel("DNA position")
    plt.ylabel("Lagging strand GC-skew : WS = 0.2 %")
    ##Fourth panel : Lagging strand cumulative GC-skew
    plt.subplot(224)
    plt.plot(range(len(c4)), c4)
    plt.xlabel("DNA position")
    plt.ylabel("WS = 0.2 % ")
    ##Save the whole figure
    plt.savefig("fig"+str(idx)+".png")
    plt.show()
    return std,mean,c2


# For selective pressure study (not applied here)
def PullCodon(skews,position):
    skewed = []
    for i in range(len(skews)):
        if skews[i]>0:
            #print ("La position du codon skewed est : ",position[i])
            skewed.append(position[i])
    print(skewed)


# For gene_position pulling (not applied here)
def ParseAnnot(line):
    gene_position = []
    for i in line:
        for j in range(len(i)):
            if i[j:j+12] =="assembly_gap":
                for m in range(j+12,j+24):
                    #print(i[m])
                    if i[m]=='.':
                        gene_position.append(i[j+12:m])

    return gene_position


# Pull random line lists from texte files
def getData_aleatoire(path):
    # Open the texte file
    genome_file = open(path, "r")
    # Get file_size
    genomeSize = os.stat(path).st_size
    # Get file lines
    lines = genome_file.readlines()
    lines_list = [line for line in lines]
    # Remove the first line
    return lines_list,genomeSize


# Pull random seqs from random line lists
def getRandom(random_seqs_paths,k):
        #### Pull random seqs :
    lines_list, genome_size = getData_aleatoire(random_seqs_paths[k])
    seqs = [lines_list[i] for i in range(len(lines_list)) if str(lines_list[i][0])!=">"]
    return seqs


# Calculate DN distribution
def DN_distribution(seqs):
    occ = [] #We create a DN occurence list
    total = [] #We create a total DN occurence list

    #We pull AT and GC content for each seqs
    for i in seqs:
        DN = {"aa":0, "at":0, "ag":0, "ac":0, "tt":0, "ta":0, "tg":0, "tc":0, "gg":0, "ga":0, "gt":0, "gc":0, "cc":0, "ca":0, "ct":0, "cg":0}
        nA = i.count("a")
        nT = i.count("t")
        nG = i.count("g")
        nC = i.count("c")
        try:
            AT_content = nA+nT/(nA+nT+nG+nC)*100
            GC_content = nG+nC/(nA+nT+nG+nC)*100
        except ZeroDivisionError:
            pass
        # We compute DN occurence for each random seq
        if len(i)>7:
            for j in range(len(i)-2):
                if i[j:j+2] in DN:
                    DN[i[j:j+2]]= DN[i[j:j+2]]+1

            occ.append(DN)
    # Then, we compute DN frequencies for each random seq
    for i in occ:
        DN = {"aa": 0, "at": 0, "ag": 0, "ac": 0, "tt": 0, "ta": 0, "tg": 0, "tc": 0, "gg": 0, "ga": 0, "gt": 0,"gc": 0, "cc": 0, "ca": 0, "ct": 0, "cg": 0}
        total.append(sum(i.values()))
    freq = [ [round(m/total[i]*100) for m in occ[i].values()] for i in range(len(occ)) ]
    return freq


# Metafunction that get, parse  calculate GC_skews and plot it giving biological sequence txt files
def analyze(Genomes,idx,label,WS):
    # Get lines_lists and genomes sizes
    lines_list,genome_size = getData(Genomes)
    # Get parsed sequences
    genome = parse(lines_list)
    # Get gene positions from txt annotations (# not applied here)
    #line2,genome_size = getData("Bacillius_anthracis\Bacilius_anthracis_annotations.txt")
    #gene_position = ParseAnnot(line2)
    # Get DNA content
    DN_content(genome,label)
    print("Genome size: ", len(genome), "bp")
    print("#####################################")
    # Plot skews
    std,mean,spectra = plot_skew(genome,label,idx,WS)
    return std,mean,genome,spectra


# Pipeling function
def GC_skew_analysis(Genomes,labels, WS,Species):

    idx = 1
    stds = []
    means = []
    UC_genomes = []
    spectrums = []
    # Plot GC_skew analysis and statistical description
    for i in range(len(Genomes)):
        std, mean, genome,spectra = analyze(Genomes[i], idx, labels[i],WS)
        UC_genomes.append("".join(genome))
        stds.append(std)
        means.append(mean)
        spectrums.append(spectra)
        idx = idx + 1
    plot_stats(stds,labels)
    return stds, means, UC_genomes, spectrums


# function to get lower cased genomes copies
def getLower(UC_genomes):
    processed_genome = []
    for i in UC_genomes:
        curr = []
        for j in range(len(i)):
            curr.append(i[j].lower())
        processed_genome.append("".join(curr))
    return processed_genome

# Plot GC skew stats
def plot_stats(stds,labels):

    plt.figure(1, figsize=(20, 16))
    plt.bar(range(len(stds)),stds,width=0.7)
    plt.ylabel('Std')
    plt.title('GC-skew std')
    plt.xticks(range(len(stds)),labels)
    plt.tick_params(axis='x', labelsize=5)
    plt.savefig("Std_measurements.png")
    plt.show()

# Pipelining function for DN distribution analysis
def Genome_signature_analysis(labels,random_seqs_paths,UC_genomes,n_random_sets):
    ### Process organisms seqs to ensure that characters are lower cased (LC needed for DN_distribution function)
    processed_genome = getLower(UC_genomes)

    ### Compute biological freqs
    biological_freqs = DN_distribution(processed_genome)

    #Set an iterator in order to walk through the biological_freqs list (ie: for each sequence)
    iter = 0

    for genome_idx in range(0,len(random_seqs_paths),n_random_sets): # for each random sets (ie: for each random species sets containing sets of sizes 10,20,40 and 80)
        replicates = [DN_distribution(getRandom(random_seqs_paths,random_replicate)) for random_replicate in range(genome_idx,genome_idx+ 4) ] #Compute DN distribution
        random_mean = []
        random_std = []
        replicates_set = [replicates[DN_frequencies][DN_frequency] for DN_frequencies in range(len(replicates)) for DN_frequency in range(len(replicates[DN_frequencies]))] #Gather frequencies

        for DN_frequency in range(len(replicates_set[0])):# for each DN frequency
            random_mean.append(np.mean([replicates_set[replicate][DN_frequency] for replicate in range(len(replicates))])) # Compute the mean
            random_std.append(np.std([replicates_set[replicate][DN_frequency] for replicate in range(len(replicates))])) # and the std

        ###Plot the figure
        fig = plt.figure(1, figsize=(20, 16))
        ax = fig.add_subplot(111)
        Labels = ("AA", "AT", "AG", "AC", "TT", "TA", "TG", "TC", "GG", "GA", "GT", "GC", "CC", "CA", "CT", "CG")
        ratios = [biological_freqs[iter][m]/random_mean[m] for m in range(len(random_mean))] #Calculate ratios between biological freq means and random freqs
        Labels = [x for _, x in sorted(zip(ratios, Labels))]
        ratios = sorted(ratios)
        barlist = plt.bar(range(len(ratios)),ratios,width=0.7)
        for i in range(len(ratios)):
            if ratios[i]>1.2:
                barlist[i].set_color('r')
            if ratios[i]<0.8:
                barlist[i].set_color('g')
        red_patch = mpatches.Patch(color='red', label='ratio > 1.2')
        blue_patch = mpatches.Patch(color='blue', label='ratio unchanged')
        green_patch = mpatches.Patch(color='green', label='ratio <0.8')
        plt.legend(handles=[red_patch,blue_patch,green_patch])
        plt.ylabel('Random vs observed dinucleotide frequencies ratios')
        plt.title(str(labels[iter])+' dinucleotide frequencies ratio')
        plt.xticks(range(len(Labels)), Labels)
        plt.savefig(str(labels[iter])+"_dinucleotide_frequencies_vs_random.png")
        plt.show()
        iter = iter+1

    ###Stastitical tests
    for i in range(len(biological_freqs)):
        pval_list=[]
        for j in range(len(biological_freqs)):
            pval = stats.kruskal(biological_freqs[i], biological_freqs[j])
            pval_list.append(pval[1])
            print(labels[i]+" vs "+labels[j]," => p-value =  ",pval[1])
        try:
            plt.ylim(0, max(pval_list))
        except ValueError:
            plt.ylim(0, 1)
        fig = plt.figure(1, figsize=(20, 16))
        ax = fig.add_subplot(111)
        barlist = plt.bar(range(len(pval_list)), pval_list,width=0.7)
        barlist[0].set_color('r')
        barlist[1].set_color('r')
        barlist[2].set_color('b')
        barlist[3].set_color('b')
        barlist[4].set_color('b')
        barlist[5].set_color('g')
        barlist[6].set_color('g')
        barlist[7].set_color('g')
        barlist[8].set_color('g')
        barlist[9].set_color('g')
        barlist[10].set_color('g')
        barlist[11].set_color('g')
        red_patch = mpatches.Patch(color='red', label='Prokaryota')
        blue_patch = mpatches.Patch(color='blue', label='Cyanobacteria')
        green_patch = mpatches.Patch(color='green', label='Eukaryota')
        plt.legend(handles=[red_patch, blue_patch, green_patch])
        plt.ylabel("p_value")
        plt.title("Statistical tests for genome signature correlation : " + labels[i])
        x = plt.xticks(range(len(biological_freqs)), labels)
        plt.tick_params(axis='x', labelsize=4)
        plt.savefig(str(labels[i]) + "Genome_signatures_statistical.png")
        plt.show()

def GC_skew_statistical_analysis(spectrums,labels):
    for i in range(len(spectrums)):
        pval_list = []
        for j in range(len(spectrums)):
            #corr1 = signal.correlate(spectrums[i], spectrums[i])
            #corr2 = signal.correlate(spectrums[i], spectrums[j])
            # Shapiro wilk normal distribution test: #
            t_shapiro1 = shapiro(spectrums[i])              #
            t_shapiro2 = shapiro(spectrums[j])              ####Functions were disabled to reduce computation time
            # Bartlett test for variance homogeneity #
            t_bartlett = bartlett(spectrums[i], spectrums[j])     #

            if t_shapiro1[1]<0.01 and t_shapiro2[1]<0.01 and t_bartlett[1]<0.01 : ### If datas follow normal distribution and have homogene variances
                pval1 = f_oneway(spectrums[i], spectrums[j]) # Do a 2 samples indep t-test
                if np.isinf(pval1[1]) == True or pval1[1] == 'Inf': # if p-value is infinite
                    pval_list.append(700) # use max value
                else:
                    pval_list.append(abs(np.log(pval1[1]))) #compute absolute log p-value
            else:
                pval1 = kruskal(spectrums[i], spectrums[j]) # Do a mann and whitney test
                if np.isinf(pval1[1]) == True or pval1[1] == 'Inf':
                    pval_list.append(700)
                else:
                    pval_list.append(abs(np.log(pval1[1])))

            print(labels[i], "vs", labels[j],pval1[1])
        try:
            plt.ylim(0, max(pval_list))
        except ValueError:
            plt.ylim(0, 700)
        fig = plt.figure(1, figsize=(20, 16))
        ax = fig.add_subplot(111)
        barlist = plt.bar(range(len(pval_list)), pval_list,width=0.7)
        barlist[0].set_color('r')
        barlist[1].set_color('r')
        barlist[2].set_color('b')
        barlist[3].set_color('b')
        barlist[4].set_color('b')
        barlist[5].set_color('g')
        barlist[6].set_color('g')
        barlist[7].set_color('g')
        barlist[8].set_color('g')
        barlist[9].set_color('g')
        barlist[10].set_color('g')
        barlist[11].set_color('g')
        red_patch = mpatches.Patch(color='red', label='Prokaryota')
        blue_patch = mpatches.Patch(color='blue', label='Cyanobacteria')
        green_patch = mpatches.Patch(color='green', label='Eucaryota')
        plt.legend(handles=[red_patch, blue_patch, green_patch])
        plt.ylabel("Absolute log p_value")
        plt.title("GC-skew profil correlation between species : " + labels[i])
        x = plt.xticks(range(len(spectrums)), labels)
        plt.tick_params(axis='x', labelsize=4)
        plt.savefig(str(labels[i]) + "GC_skew_statistical.png")
        plt.show()