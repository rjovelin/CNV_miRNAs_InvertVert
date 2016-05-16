from CNV_miRNAs_TargetScan import *
from scipy import stats



def mirna_regulation_UTR_length_controlled(CNV_miRNA_file, UTR_file, species, regulator):
    '''
    (file, file, str, str) -> tuple of tuples
    Perform a Wilcoxon rank sum test for regulator mirnas or sites normalized by the length of the UTR
    and return a tuple with the sample size, mean and sem of the regulator variable for CNV genes, a tuple with the
    sample size, the mean and sem for non_CNV genes and a tuple with the Wilcoxon z-score and its p-value
    '''

    # open file for reading
    CNVs = open(CNV_miRNA_file, 'r')
    # skip header
    CNVs.readline()

    # make a dictionnary with the CNV status for transcript {transcript: [gene, N_mirnas, N_sites, CNV_status]}
    target_transcripts = {}

    # go through the file and populate the dict
    for line in CNVs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            transcript = line[1]
            mirnas = int(line[2])
            sites = int(line[3])
            CNV_status = line[-1]
            target_transcripts[transcript] = [gene, mirnas, sites, CNV_status]

    # close file after reading
    CNVs.close()

    # get the sequence of all UTRs for the given species {transcript : UTR}
    all_UTR = UTR_sequence(UTR_file, species)


    # make a list of mirna_regulation for CNV genes and one for non-CNV genes
    CNV_regulation = []
    nonCNV_regulation = []

    # populate the lists with the number of mirnas (or sites) / length of UTR
    for transcript in target_transcripts:
        UTR_length = len(all_UTR[transcript].replace('-', ''))
        if regulator == 'mirnas':
            N_regulators = target_transcripts[transcript][1]
        elif regulator == 'sites':
            N_regulators = target_transcripts[transcript][2]
        if target_transcripts[transcript][-1] == 'CNV':
            CNV_regulation.append(N_regulators / UTR_length)
        elif target_transcripts[transcript][-1] == 'non-CNV':
            nonCNV_regulation.append(N_regulators / UTR_length)


    # perform a Wilcoxon rank sum test
    diff = stats.ranksums(CNV_regulation, nonCNV_regulation)
    diff1 = stats.ttest_ind(CNV_regulation, nonCNV_regulation)
    diff2 = stats.ks_2samp(CNV_regulation, nonCNV_regulation)

    # compute the mean and sem for CNV and non-CNV genes
    mean_CNV = stats.tmean(CNV_regulation)
    sem_CNV = stats.sem(CNV_regulation)
    mean_nonCNV = stats.tmean(nonCNV_regulation)
    sem_nonCNV = stats.sem(nonCNV_regulation)

    
    return (len(CNV_regulation), mean_CNV, sem_CNV), (len(nonCNV_regulation), mean_nonCNV, sem_nonCNV), diff, diff1, diff2




##def sort_UTR_length_CNV_status(CNV_miRNA_file, UTR_file, species, regulator):
##
##    ######### this function is deprecated #################
##
##
##    
##
##    # open file for reading
##    CNVs = open(CNV_miRNA_file, 'r')
##    # skip header
##    CNVs.readline()
##
##    # make a dictionnary with the CNV status for transcript {transcript: [gene, N_mirnas, N_sites, CNV_status]}
##    target_transcripts = {}
##
##    # go through the file and populate the dict
##    for line in CNVs:
##        line = line.rstrip()
##        if line != '':
##            line = line.split()
##            gene = line[0]
##            transcript = line[1]
##            mirnas = int(line[2])
##            sites = int(line[3])
##            CNV_status = line[-1]
##            target_transcripts[transcript] = [gene, mirnas, sites, CNV_status]
##
##    # close file after reading
##    CNVs.close()
##
##    # get the sequence of all UTRs for the given species {transcript : UTR}
##    all_UTR = UTR_sequence(UTR_file, species)
##
##    # make a dict with the UTR length of the target transcripts
##    targets_UTR = {}
##    for transcript in target_transcripts:
##        targets_UTR[transcript] = len(all_UTR[transcript])
##
##    # make a list of UTR lengths
##    UTR_lengths = []
##    for transcript in targets_UTR:
##        UTR_lengths.append(targets_UTR[transcript])
##
##    # get the maximum UTR length
##    longest = max(UTR_lengths)
##
##    # make a dict to hold the transctip names that have UTR within a 100 bp range
##    # use a index  = length(UTR) // 100 to find which transcripts should be grouped together
##
##    transcripts_length = {}
##    for i in range(0, longest, 100):
##        transcripts_length[int(i // 100)] = []
##
##    # populate the dict with the name of the transcripts
##    for transcript in target_transcripts:
##        position = int(targets_UTR[transcript]) // 100
##        transcripts_length[position].append(transcript)
##
##    # create a dict to hold the number of mirnas or sites for CNV genes
##    CNV_length = {}
##    # use the same keys as the dictionnary holding the transcript names for UTR length range
##    for i in transcripts_length:
##        CNV_length[i] = 0
##
##    # create a dict to hold the number of mirnas or sites for non-CNV genes
##    nonCNV_length = {}
##    # use the same keys as the dictionnary holding the transcript names for UTR length range
##    for i in transcripts_length:
##        nonCNV_length[i] = 0
##
##
##    # get the number of mirnas or sites for each UTR length range for CNV genes and non-CNV genes
##    for i in transcripts_length:
##        # count mirnas or sites per gene
##        if regulator == 'mirnas':
##            # reinitialise the counters
##            N_mirnas_CNV = 0
##            N_mirnas_nonCNV = 0
##            # reinitialize the transcript counter
##            N_CNVs = 0
##            N_nonCNVs = 0
##            # check that transcripts are within that UTR length range
##            if len(transcripts_length[i]) != 0:
##                # count the mean number of mirnas per gene
##                for transcript in transcripts_length[i]:
##                    if target_transcripts[transcript][-1] == 'CNV':
##                        N_mirnas_CNV += target_transcripts[transcript][1]
##                        N_CNVs +=1
##                    elif target_transcripts[transcript][-1] == 'non-CNV':
##                        N_mirnas_nonCNV += target_transcripts[transcript][1]
##                        N_nonCNVs += 1
##                # populate the dicts for CNVs and non-CNVs with the average number of mirnas
##                if N_CNVs != 0:
##                    CNV_length[i] = N_mirnas_CNV / N_CNVs
##                if N_nonCNVs != 0:
##                    nonCNV_length[i] = N_mirnas_nonCNV / N_nonCNVs
##        elif regulator == 'sites':
##            # reinitialise the counters
##            N_sites_CNV = 0
##            N_sites_nonCNV = 0
##            # reinitialize the transcript counter
##            N_CNVs = 0
##            N_nonCNVs = 0
##            # check that transcripts are within that UTR length range
##            if len(transcripts_length[i]) != 0:
##                # count the mean number of sites per gene
##                for transcript in transcripts_length[i]:
##                    if target_transcripts[transcript][-1] == 'CNV':
##                        N_sites_CNV += target_transcripts[transcript][2]
##                        N_CNVs +=1
##                    elif target_transcripts[transcript][-1] == 'non-CNV':
##                        N_sites_nonCNV += target_transcripts[transcript][2]
##                        N_nonCNVs += 1
##                # populate the dicts for CNVs and non-CNVs with the average number of mirnas
##                if N_CNVs != 0:
##                    CNV_length[i] = N_sites_CNV / N_CNVs
##                if N_nonCNVs != 0:
##                    nonCNV_length[i] = N_sites_nonCNV / N_nonCNVs
##                        
##
##    # make a list of index
##    lower_bound_CNVs = [i for i in CNV_length]
##    lower_bound_nonCNVs = [i for i in nonCNV_length]
##    assert len(lower_bound_CNVs) == len(lower_bound_nonCNVs)
##
##    # sort lists
##    lower_bound_CNVs.sort()
##    lower_bound_nonCNVs.sort()
##
##    # make list with mirnas or sites counts ordered by index of UTR length range
##    CNV_regulation = []
##    for i in lower_bound_CNVs:
##        CNV_regulation.append(CNV_length[i])
##    nonCNV_regulation = []
##    for i in lower_bound_nonCNVs:
##        nonCNV_regulation.append(nonCNV_length[i])
##
##    # perform the Kolmogorov-2 sample test between the distributions of CNV and non-CNV regulation
##    KS_2S_test = stats.ks_2samp(CNV_regulation, nonCNV_regulation)
##
##    return CNV_regulation, nonCNV_regulation, KS_2S_test





    









    









    
