from CNV_miRNAs_TargetScan import *
from cnv_mirnas_miranda import *


def get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file):
    '''
    (file) -> tuple
    Returns a 6-item tuple containing the lists of the number of miRNA regulators,
    the number of miRNA binding sites, and the number of sites per miRN for CNV and non-CNV genes
    '''

    # open CNV_miRNA file
    CNV_miRNA = open(CNV_miRNA_file, 'r')

    # create lists to store the values of the different metrics for CNV and non-CNV genes
    CNV_N_mirnas = []
    CNV_N_sites = []
    CNV_ratio = []
    nonCNV_N_mirnas = []
    nonCNV_N_sites = []
    nonCNV_ratio = []

    # read file and store the the values in the lists
    CNV_miRNA.readline()
    for line in CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[-1] == 'CNV':
                CNV_N_mirnas.append(int(line[2]))
                CNV_N_sites.append(int(line[3]))
                CNV_ratio.append(float(line[4]))
            elif line[-1] == 'non-CNV':
                nonCNV_N_mirnas.append(int(line[2]))
                nonCNV_N_sites.append(int(line[3]))
                nonCNV_ratio.append(float(line[4]))
        
    CNV_miRNA.close()
    return CNV_N_mirnas, nonCNV_N_mirnas, CNV_N_sites, nonCNV_N_sites, CNV_ratio, nonCNV_ratio


def test_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file):
    '''
    (file) -> tuple
    Performs a Wilcoxon rank sum test between CNV and non-CNV genes of
    the number of miRNA regulators, the number of miRNA binding sites, and the number of sites per miRNA
    Return a tuple of 2-item tuple containing the z-value and the p-value
    '''

    # get the number of miRNAs, sites and sites per miRNA for CNV and non-CNV genes
    miRNA_regulation = get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
    CNV_N_mirnas = miRNA_regulation[0]
    nonCNV_N_mirnas = miRNA_regulation[1]
    CNV_N_sites = miRNA_regulation[2]
    nonCNV_N_sites = miRNA_regulation[3]
    CNV_ratio = miRNA_regulation[4]
    nonCNV_ratio = miRNA_regulation[5]

    # compute the Wilcoxon rank sum tests
    from scipy import stats
    diff_N_mirnas = stats.ranksums(CNV_N_mirnas, nonCNV_N_mirnas)
    diff_N_sites = stats.ranksums(CNV_N_sites, nonCNV_N_sites)
    diff_ratio = stats.ranksums(CNV_ratio, nonCNV_ratio)

    return diff_N_mirnas, diff_N_sites, diff_ratio


def get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file):
    '''
    (file, str, file) -> tuple
    Returns a 2-item tuple with the lists of UTR length for CNV and non-CNV genes recorded in CNV_miRNA_file
    '''

    # get the sequence of all UTRs for the given species
    all_UTR = UTR_sequence(UTR_file, species)

    # make a dictionnary of CNV status for each transcript in the CNV_miRNA_file
    UTR_target_CNV_status = {}
    CNV_miRNA = open(CNV_miRNA_file, 'r')
    CNV_miRNA.readline()
    for line in CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            UTR_target_CNV_status[line[1]] = line[-1]

    # make lists to the store the UTR length of CNV and non-CNV target transcripts
    CNV_UTR = []
    nonCNV_UTR = []

    # partition the UTR length of transcripts target based on CNV status
    for UTR in all_UTR:
        if UTR in UTR_target_CNV_status:
            if UTR_target_CNV_status[UTR] == 'CNV':
                CNV_UTR.append(len(all_UTR[UTR]))
            elif UTR_target_CNV_status[UTR] == 'non-CNV':
                nonCNV_UTR.append(len(all_UTR[UTR]))

    CNV_miRNA.close()
    return CNV_UTR, nonCNV_UTR
       
    
def test_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file):
    '''
    (file, str, file) -> tuple
    Performs a Wilcoxon rank sum test of the UTR length difference between CNV and non-CNV genes recorded in CNV_miRNA_file for a given species
    Returns a tuple with the z-value and the p-value
    '''

    # get the length of UTR for CNV and non-CNV target genes
    UTR_length = get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)
    CNV_UTR = UTR_length[0]
    nonCNV_UTR = UTR_length[1]
    
    # performs the Wilcoxon rank sum test
    from scipy import stats
    diff_UTR_length = stats.ranksums(CNV_UTR, nonCNV_UTR)

    return diff_UTR_length


def compute_mean_std_error(L):
    '''
    (list) -> tuple
    Returns a tuple containing the mean and the standard error of a collection of values in the list L
    Pre-condition: the values in L are floats and/or integers
    '''

    # verify the pre-condition
    for item in L:
        try:
            item + 1
        except:
            print('values in L need to be intergers and/or floats')

    import math

    # compute the mean
    total = 0
    for item in L:
        total += item
    mean = total/ len(L)

    # compute the stand error of the mean
    total_diff = 0
    for item in L:
        total_diff += (item - mean)**2
    std_dev = math.sqrt(total_diff / len(L))
    std_error = std_dev / math.sqrt(len(L))

    return (mean, std_error)



def get_UTR_length_miRNA_regulation(UTR_file, species, CNV_miRNA_file):
    '''
    (file, str, file) -> list
    Returns a list of lists containing the length of target UTRs from CNV_miRNA_file for a given species with the same order
    as the lists of the number of miRNAs, binding sites and ratio extracted from the same file
    '''

    # get the sequence of all UTRs for the given species
    all_UTR = UTR_sequence(UTR_file, species)

    # make a dictionnary with information for each target transcript
    UTR_target = {}
    CNV_miRNA = open(CNV_miRNA_file, 'r')
    CNV_miRNA.readline()
    for line in CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            UTR_target[line[1]] = [int(line[2]), int(line[3]), float(line[4])]

    # make a list to the store the UTR length of miRNA targets 
    UTR_length = []
    
    # create lists to store the number of miRNAs, binding sites and sites per miRNA for the targets
    N_mirnas = []
    N_sites = []
    ratio = []
    
    # populate the different lists keeping the same order within each list
    for transcript in UTR_target:
            UTR_length.append(len(all_UTR[transcript]))
            N_mirnas.append(UTR_target[transcript][0])
            N_sites.append(UTR_target[transcript][1])
            ratio.append(UTR_target[transcript][2])
                    
    CNV_miRNA.close()
    return [UTR_length, N_mirnas, N_sites, ratio]


def save_UTR_mirnas_to_file(UTR_file, species, CNV_miRNA_file, outputfile):
    '''
    (file, str, file, file) -> file
    Gets the list of UTR length and the number of mirna regulators for each transcript in CNV_mirna_file and saves them in a file for plotting the correlation
    between the number of mirna regulators and the length of the UTR
    '''

    # get the lists of UTR length and number of mirnas
    utr_mirna = get_UTR_length_miRNA_regulation(UTR_file, species, CNV_miRNA_file)
    UTR_length = utr_mirna[0]
    N_mirnas = utr_mirna[1]
    
    # write the content of each list to file
    newfile = open(outputfile, 'w')
    newfile.write('N_mirnas' + '\t' + 'UTR_length' + '\n')
    for i in range(len(UTR_length)):
        newfile.write(str(N_mirnas[i]) + '\t' + str(UTR_length[i]) + '\n')

    newfile.close()


def correlation_UTR_length_miRNA_regulation(UTR_file, species, CNV_miRNA_file):
    '''
    (file, str, file) -> tuple
    Performs a Spearman's rank correlation between the UTR length of each target and:
    number the miRNAs, number of sites and ratio of sites per miRNA of the same target
    Returns a tuple of tuple with each inner tuple containing the rho and p values
    '''

    from scipy import stats

    # get the list of UTR length, and N mirnas, sites and ratio
    regulated_targets = get_UTR_length_miRNA_regulation(UTR_file, species, CNV_miRNA_file)

    # unpack the list of lists to get the lists of values for each variable
    UTR_length, N_mirnas, N_sites, ratio = regulated_targets

    # perform Sperman's rank correlation between UTR length and number of miRNAs
    rho_mirnas, p_mirnas = stats.spearmanr(UTR_length, N_mirnas)

    # perform Sperman's rank correlation between UTR length and number of sites
    rho_sites, p_sites = stats.spearmanr(UTR_length, N_sites)

    # perform Sperman's rank correlation between UTR length and ratio of sites per mirna
    rho_ratio, p_ratio = stats.spearmanr(UTR_length, ratio)

    return (rho_mirnas, p_mirnas), (rho_sites, p_sites), (rho_ratio, p_ratio)


def print_results_mean_comparison_CNV_nonCNV(CNV_miRNA_file, file_type, **species_UTR):
    '''
    (file, str) -> None
    Print the results of

    '''

    if file_type == 'miranda' or file_type == 'targetscan':
        test_mirna = test_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        list_mirnas = get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        
        wilcoxon = ['mirnas_cnv_noncnv:', 'sites_cnv_noncnv:', 'ratio_cnv_noncnv:']
        print('Wilcoxon rank sum tests:')
                
        for i in range(len(test_mirna)):
            print(wilcoxon[i], 'z-score = ', test_mirna[i][0], 'p = ', test_mirna[i][1], sep = '\t')
            
        print('\n')
                
        print('N CNV genes:' + '\t' + str(len(list_mirnas[0])))
        print('N non-CNV genes:' + '\t' + str(len(list_mirnas[1])))

        print('\n')
        
        print('\t' + 'mean' + '\t' + 'std_error')
        headers = ['CNV_N_mirnas:', 'non-CNV_N_mirnas:',
                   'CNV_N_sites:', 'non-CNV_N_sites:',
                   'CNV_sites_per_mirna:', 'non-CNV_sites_per_mirna:']
        
        for i in range(len(list_mirnas)):
            mean_stderr = compute_mean_std_error(list_mirnas[i])
            print(headers[i], mean_stderr[0], mean_stderr[1], sep = '\t')

        print('\n')
        
    if file_type == 'targetscan':
        file_kw = sorted(species_UTR)
        UTR_file = species_UTR[file_kw[0]]
        species = species_UTR[file_kw[1]]
        utr_test = test_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)
        print('Wilcoxon rank sum tests:')
        print('test UTR length:', 'z-score:', utr_test[0], 'p', utr_test[1], sep = '\t')

        print('\n')

        utr_length = get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)

        print('\t' + 'mean' + '\t' + 'std_error')
        utr_header = ['CNV_N_utr:', 'non-CNV_utr:']
        
        for i in range(len(utr_length)):
            mean_stderr = compute_mean_std_error(utr_length[i])
            print(utr_header[i], mean_stderr[0], mean_stderr[1], sep = '\t')


def human_single_study(CNV_file, all_gene_file, UTR_file, summary_counts, species, conservation_sites, site_score, bonferroni,
                       file_type, miRNA_family_info, CNV_miRNA_file):
    '''
    (file, file, file, file, str, str, str, str, str, str, file, file) -> dict
    Test for each individual study differences in mean number of mirnas, sites and sites per mirna and UTR length between CNV genes
    and non-CNV genes. Returns a dictionnary with study as key and a list of outcomes as value (significantly lower, greater, not_diff)     
    '''

    
    # get the set of studies
    studies = set()
    cnv = open(CNV_file, 'r')
    cnv.readline()
    for line in cnv:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            studies.add(line[6])
    cnv.close()

    # create dictionnary to store the studies corresponding to the different outcome
    sorted_studies = {}
    
    # get the gene ID for each target transcript
    transcripts = transcript_gene_pairs(summary_counts)
    
    # get the CNV status and sort the valid genes
    valid_genes = sort_valid_human_genes(all_gene_file)
    
    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

    # for each study
    for study in studies:
        # get the set of cnv genes
        CNV_genes = human_CNV_genes_single_study(CNV_file, all_gene_file, study)
        print(len(CNV_genes), study, sep = '\t')
        # create a dictionnary to store info about miRNA regultors for each gene if file == targescan
        gene_regulated = {}
        # make a dictionnary {gene1:[N_mirnas, N_sites, ratio, transcript, CNV] if file == miranda
        regulated_genes = {}

        if file_type == 'targetscan':
            # get the mirna information
            regulators = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)

            # write the mirna info to file; using a single transcript per gene
            # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, CNV_status]}
            
            for transcript_ID in regulators:
                gene = transcripts[transcript_ID]
                if gene not in gene_regulated:
                    gene_regulated[gene] = regulators[transcript_ID] # get the info for a single transcript
                    gene_regulated[gene].append(transcript_ID) # add the transcript ID

            # add the CNV status
            for gene in gene_regulated:
                if gene in CNV_genes:
                    gene_regulated[gene].append('CNV')
                else:
                    gene_regulated[gene].append('non-CNV')

            # remove non-valid target genes
            to_remove = []
            for gene in gene_regulated:
                if gene not in valid_genes:
                    to_remove.append(gene)
            for gene in to_remove:
                del gene_regulated[gene]

            # write to file if > 500 CNV targets
            if len(CNV_genes) > 500:
                newfile = open(CNV_miRNA_file, 'w')
                newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
                for gene in gene_regulated:
                    newfile.write(gene + '\t')
                    newfile.write(gene_regulated[gene][-2] + '\t')
                    N_mirnas = 0
                    N_sites = 0
                    for pair in gene_regulated[gene][:-2]:
                        N_mirnas += len(family[pair[0]])
                        N_sites += pair[1]
                    newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/ N_mirnas) + '\t' + gene_regulated[gene][-1] + '\n')

                newfile.close()
                print('file written for' + '\t' + study)

        elif file_type == 'miranda':

            # get the transcript : gene pairs
            transcripts_genes = transcripts_gene_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)
    
            # get the number of regulators for each transcript
            regulated_transcripts = miRNA_regulation_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)

            # populate dictionary regulated_genes
            for transcript in regulated_transcripts:
                gene = transcripts_genes[transcript]
                if gene not in regulated_genes:
                    N_mirnas = len(regulated_transcripts[transcript])
                    N_sites = 0
                    for mirna in regulated_transcripts[transcript]:
                        N_sites += regulated_transcripts[transcript][mirna]
                    regulated_genes[gene] = [N_mirnas, N_sites, N_sites/N_mirnas]
                    regulated_genes[gene].append(transcript)

            # add the CNV status to each gene in regulated_genes
            for gene in regulated_genes:
                if gene in valid_genes:
                    if gene in CNV_genes:
                        regulated_genes[gene].append('CNV')
                    else:
                        regulated_genes[gene].append('non-CNV')
    
            # remove target genes that are not valid
            to_remove = []
            for gene in regulated_genes:
                if gene not in valid_genes:
                    to_remove.append(gene)
            if len(to_remove) != '':
                for gene in to_remove:
                    del regulated_genes[gene]
    
            # write to file if > 500 CNV targets
            if len(CNV_genes) > 500:
                newfile = open(CNV_miRNA_file, 'w')
                newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
                for gene in regulated_genes:
                    newfile.write(gene + '\t' + regulated_genes[gene][-2] + '\t' + str(regulated_genes[gene][0]) + '\t' + str(regulated_genes[gene][1]) + '\t'
                                  + str(regulated_genes[gene][2]) + '\t' + regulated_genes[gene][-1] + '\n')

                newfile.close()
                print('file written for' + '\t' + study)

        # grab the content of the CNV_miRNA_file if > 500 CNV targets
        if len(CNV_genes) > 500:
            list_mirnas = get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        
            # get the lists of UTR length
            CNV_UTR, nonCNV_UTR = get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)

            # test that UTR length of CNV and nonCNV genes is different, keep the p-value
            p_UTR = test_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)[1]

            # test that mean number of mirnas, sites and sites per mirna are different netween CNV and nonCNV genes
            test_mirna = test_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
            p_mirna, p_sites, p_ratio = test_mirna[0][1], test_mirna[1][1], test_mirna[2][1]

            # compute means for the CNV and non-CNV genes for the 3 variables and for UTR length
            CNV_mean_mirnas = compute_mean_std_error(list_mirnas[0])[0]
            nonCNV_mean_mirnas = compute_mean_std_error(list_mirnas[1])[0]
            CNV_mean_sites = compute_mean_std_error(list_mirnas[2])[0]
            nonCNV_mean_sites = compute_mean_std_error(list_mirnas[3])[0]
            CNV_mean_ratio = compute_mean_std_error(list_mirnas[4])[0]
            nonCNV_mean_ratio = compute_mean_std_error(list_mirnas[5])[0]
            CNV_mean_UTR = compute_mean_std_error(CNV_UTR)[0]
            nonCNV_mean_UTR = compute_mean_std_error(nonCNV_UTR)[0]

            # ask if the test is significant and which of CNV of non-CNV genes has greater mean
            if bonferroni == 'yes':
                threshold = 0.05 / 31 # 31 studies have > 500 CNV miRNA target
            elif bonferroni == 'no':
                threshold = 0.05

            # populate the dictionary
            sorted_studies[study] = []
            
            # compare mean number of mirnas
            if p_mirna < threshold:
                if CNV_mean_mirnas < nonCNV_mean_mirnas:
                    sorted_studies[study].append('mirnas_CNV_lower')
                elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                    sorted_studies[study].append('mirnas_CNV_greater')
            else:
                sorted_studies[study].append('mirnas_CNV_notdiff')
            # compare mean number of sites
            if p_sites < threshold:
                if CNV_mean_sites < nonCNV_mean_sites:
                    sorted_studies[study].append('sites_CNV_lower')
                elif CNV_mean_sites > nonCNV_mean_sites:
                    sorted_studies[study].append('sites_CNV_greater')
            else:
                sorted_studies[study].append('sites_CNV_notdiff')
            # compare mean number of sites per mirna
            if p_ratio < threshold:
                if CNV_mean_ratio < nonCNV_mean_ratio:
                    sorted_studies[study].append('ratio_CNV_lower')
                elif CNV_mean_ratio > nonCNV_mean_ratio:
                    sorted_studies[study].append('ratio_CNV_greater')
            else:
                sorted_studies[study].append('ratio_CNV_notdiff')

            # compare mean UTR length
            if p_UTR < threshold:
                if CNV_mean_UTR < nonCNV_mean_UTR:
                    sorted_studies[study].append('UTR_CNV_lower')
                elif CNV_mean_UTR > nonCNV_mean_UTR:
                    sorted_studies[study].append('UTR_CNV_greater')
            else:
                sorted_studies[study].append('UTR_CNV_notdiff')

            print('length sorted_studies' + '\t' + str(len(sorted_studies)))

    return sorted_studies


def conditional_probabilities(replicate_tests, variable):
    '''
    (dict) -> list
    Take a dictionnary with studies or replicates as key and a list of outcomes for tests of mean difference for a given variable (number of mirnas,
    sites, sites per mirna), for UTR length between CNV and non-CNV genes and returns the a list with the number of studies/replicates in which
    CNV > non-CNV for variable, for UTR and for both variable and UTR, CNV < non-CNV for variable, for UTR and for both.
    Use these numbers to compute the conditional probabilities that a larger (or lower) mean for CNV is due to a longer (or shorter) UTR
    '''

    # set the variables to compute
    regulator_CNV_greater = 0
    UTR_CNV_greater = 0
    regulator_UTR_CNV_greater = 0
    regulator_CNV_lower = 0
    UTR_CNV_lower = 0
    regulator_UTR_CNV_lower = 0

    # go through the dictionnary and update the UTR variables
    for replicate in replicate_tests:
        if replicate_tests[replicate][3] == 'UTR_CNV_greater':
            UTR_CNV_greater += 1
        elif replicate_tests[replicate][3] == 'UTR_CNV_lower':
            UTR_CNV_lower += 1
    
    # go through the dictionnary and update the values of the variables
    if variable == 'mirnas':
        for replicate in replicate_tests:
            if replicate_tests[replicate][0] == 'mirnas_CNV_greater':
                regulator_CNV_greater += 1
            elif replicate_tests[replicate][0] == 'mirnas_CNV_lower':
                regulator_CNV_lower += 1
            if replicate_tests[replicate][0] == 'mirnas_CNV_greater' and replicate_tests[replicate][3] == 'UTR_CNV_greater':
                regulator_UTR_CNV_greater += 1
            elif replicate_tests[replicate][0] == 'mirnas_CNV_lower' and replicate_tests[replicate][3] == 'UTR_CNV_lower':
                regulator_UTR_CNV_lower += 1


    elif variable == 'sites':
        for replicate in replicate_tests:
            if replicate_tests[replicate][1] == 'sites_CNV_greater':
                regulator_CNV_greater += 1
            elif replicate_tests[replicate][1] == 'sites_CNV_lower':
                regulator_CNV_lower += 1
            if replicate_tests[replicate][1] == 'sites_CNV_greater' and replicate_tests[replicate][3] == 'UTR_CNV_greater':
                regulator_UTR_CNV_greater += 1
            elif replicate_tests[replicate][1] == 'sites_CNV_lower' and replicate_tests[replicate][3] == 'UTR_CNV_lower':
                regulator_UTR_CNV_lower += 1

    elif variable == 'ratio':
        for replicate in replicate_tests:
            if replicate_tests[replicate][2] == 'ratio_CNV_greater':
                regulator_CNV_greater += 1
            elif replicate_tests[replicate][2] == 'ratio_CNV_lower':
                regulator_CNV_lower += 1
            if replicate_tests[replicate][2] == 'ratio_CNV_greater' and replicate_tests[replicate][3] == 'UTR_CNV_greater':
                regulator_UTR_CNV_greater += 1
            elif replicate_tests[replicate][2] == 'ratio_CNV_lower' and replicate_tests[replicate][3] == 'UTR_CNV_lower':
                regulator_UTR_CNV_lower += 1
    

    return [regulator_CNV_greater, UTR_CNV_greater, regulator_UTR_CNV_greater, regulator_CNV_lower, UTR_CNV_lower, regulator_UTR_CNV_lower]


def resampling_random_genes(summary_counts, species, conservation_sites, miRNA_family_info, UTR_file, all_gene_file,
                            CNV_file, replicates, ensembl_transcript_coordinates = 'zebrafish_ensembl_genes_zv8',
                            fly_gene_coordinates = 'dmel-all-r5.1.gff', zebrafish_CNV_genes_file = 'Zebrafish_CNV_genes.txt',
                            ensembl_gene_transcripts = 'Zebrafish_Gene_Transcripts_ID', CNV_miRNA_file = 'temp_outputfile.txt'):
    '''
    (file, str, str, file, file, file, file, str, int, file, file, file, file, file) -> dict, dict
    Test for each replicate differences in mean number of mirnas, sites and sites per mirna and UTR length between CNV genes
    and non-CNV genes. Returns 2 dictionnaries with study as key and a list of outcomes as value (significantly lower, greater, not_diff), one without correcting
    for multiple testing and one with Bonferroni correction
    '''

    import random
    
    # create dictionnary to store the outcomes of the tests for each replicate
    samples = {}
    corrected_samples = {}

    # get the miRNA info for each transcript
    if species == 'worm':
        regulated_transcripts = worm_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)
    elif species == 'human':
        regulated_transcripts = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)
    elif species == 'zebrafish':
        regulated_transcripts = zebrafish_miRNA_regulation_TargetScan_sites(summary_counts, species, miRNA_family_info)
    elif species == 'fly':
        regulated_transcripts = fly_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)

    # get the UTR sequences of all transcripts for each transcript
    UTR_sequences = UTR_sequence(UTR_file, species)
    
    # get the gene ID for each target transcript
    if species == 'fly':
        transcripts = fly_transcript_gene_pairs(UTR_file)
    elif species == 'worm' or species == 'human':
        transcripts = transcript_gene_pairs(summary_counts)
    elif species == 'zebrafish':
        transcripts = {}
        for transcript_ID in regulated_transcripts:
            gene_ID = transcript_ID[:transcript_ID.index('.')] # the transcript ID is not a valid ensembl ID but the gene_ID + .number
            transcripts[transcript_ID] = gene_ID      

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

    # make a reverse dictionnary with gene: [list of transcript]
    target_genes = {}
    for transcript_ID in transcripts:
        gene_ID = transcripts[transcript_ID]
        if gene_ID in target_genes:
            target_genes[gene_ID].append(transcript_ID)
        else:
            target_genes[gene_ID] = [transcript_ID]
        
    # get the set of valid genes
    if species == 'worm':
        valid_genes = sort_valid_worm_genes(all_gene_file)
    elif species == 'human':
        valid_genes = sort_valid_human_genes(all_gene_file)
    elif species == 'zebrafish':
        valid_genes = sort_valid_zebrafish_genes(ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only = 'yes')
    elif species == 'fly':
        valid_genes = sort_valid_fly_genes(fly_gene_coordinates)

    # remove genes that are not in the reference genome
    to_remove = []
    for gene_ID in target_genes:
        if gene_ID not in valid_genes:
            to_remove.append(gene_ID)
    if len(to_remove) != 0:
        for gene_ID in to_remove:
            del target_genes[gene_ID]

    # get the CNV genes
    if species == 'worm':
        CNV_genes = worm_CNV_genes(CNV_file, all_gene_file)
    elif species == 'human':
        CNV_genes = human_CNV_genes(CNV_file, all_gene_file)
    elif species == 'zebrafish':
        CNV_file = zebrafish_CNV_genes_file
        cnv = open(CNV_file, 'r')
        CNV_genes = set()
        for line in cnv:
            line = line.rstrip()
            if line != '':
                CNV_genes.add(line)
        cnv.close()
    elif species == 'fly':
        CNV_genes = set()
        cnv = open(CNV_file, 'r')
        cnv.readline()
        for line in cnv:
            line = line.rstrip()
            if line != '':
                line = line.split()
                CNV_genes.add(line[0])
                CNV_genes.add(line[1])
        cnv.close()
        
    # make separate directionnaries for targets: CNV and non-CNV
    CNV_targets = {}
    nonCNV_targets = {}

    for gene in target_genes:
        if gene in CNV_genes:
            CNV_targets[gene] = target_genes[gene]
        elif gene not in CNV_genes:
            nonCNV_targets[gene] = target_genes[gene]

    # make a list of genes for each target gene dictionnaries
    draw_from_CNV_targets = [gene for gene in CNV_targets]
    draw_from_nonCNV_targets = [gene for gene in nonCNV_targets]
    print(len(draw_from_CNV_targets))
    print(len(draw_from_nonCNV_targets))

    

    # repeat the tests for number of replicates, and store the outcomes in the dictionnary of replicates outcomes
    j = replicates
    while j != 0:
        # create a dictionnary to store the miRNA content of each transcript (1 transcript per gene) for CNV genes and for non-CNV genes
        mirna_target_CNV = {}
        mirna_target_nonCNV = {}
        
        # create a dictionnary to store mirna information for all 500 CNV genes and all 500 non-CNV genes
        regulators = {}
        
        # create a set of genes that are already picked
        already_picked = set()
        
        # pick 500 CNV genes at random
        i = 500
        while i != 0:
            # pick a random number and get the gene from the list of cnv_targets using the random index
            position = random.randint(0, len(draw_from_CNV_targets)-1)
            target = draw_from_CNV_targets[position]
            # if the gene has multiple transcripts, pick a random transcript
            if len(CNV_targets[target]) > 1:
                position = random.randint(0, len(CNV_targets[target]) -1)
                target_transcript = CNV_targets[target][position]
            else:
                target_transcript = CNV_targets[target][0]
            
            # verify that the genes were not already picked
            if target not in already_picked:
                # add the genes to set of already picked genes
                already_picked.add(target)
                #and that the transcript is in regulated_transcripts (ie. celegans includes regulation for miR* that were removed)
                if target_transcript in regulated_transcripts:
                    # use the dictionnary of transcript:[mirna info] to extract the mirna info
                    #WARNING: use list to grab the list so that regulated_transcripts is not modified
                    mirna_target_CNV[target_transcript] = list(regulated_transcripts[target_transcript]) 
                    mirna_target_CNV[target_transcript].append(target_transcript)
                    mirna_target_CNV[target_transcript].append('CNV')
                    # change i value
                    i -= 1
                
        # pick 500 non-CNV genes at random
        k = 500
        while k != 0:
            # pick a random number and get the gene from the list of non-cnv_targets using the random index
            position = random.randint(0, len(draw_from_nonCNV_targets)-1)
            target = draw_from_nonCNV_targets[position]
            # if the gene has multiple transcripts, pick a random transcript
            if len(nonCNV_targets[target]) > 1:
                position = random.randint(0, len(nonCNV_targets[target]) - 1)
                target_transcript = nonCNV_targets[target][position]
            else:
                target_transcript = nonCNV_targets[target][0]

            # verify that the genes were not already picked
            if target not in already_picked:
                # add the genes to set of already picked genes
                already_picked.add(target)
                #and that the transcript is in regulated_transcripts (ie. celegans includes regulation for miR* that were removed)
                if target_transcript in regulated_transcripts:
                    # use the dictionnary of transcript:[mirna info] to extract the mirna info
                    #WARNING: use list to grab the list so that regulated_transcripts is not modified
                    mirna_target_nonCNV[target_transcript] = list(regulated_transcripts[target_transcript])
                    mirna_target_nonCNV[target_transcript].append(target_transcript)
                    mirna_target_nonCNV[target_transcript].append('non-CNV')
                    # change k value
                    k -= 1
        
        # populate the dictionnary regulators with mirna information, transcript_ID and CNV status
        # write content to temporary file       
        for target_transcript in mirna_target_CNV:
            if species == 'worm' or species == 'human' or species == 'fly':
                gene = transcripts[target_transcript]
            elif species == 'zebrafish':
                gene = target_transcript[:target_transcript.index('.')]
            regulators[gene] = mirna_target_CNV[target_transcript]
        for target_transcript in mirna_target_nonCNV:
            if species == 'worm' or species == 'human' or species == 'fly':
                gene = transcripts[target_transcript]
            elif species == 'zebrafish':
                gene = target_transcript[:target_transcript.index('.')]
            regulators[gene] = mirna_target_nonCNV[target_transcript]
       
        newfile = open(CNV_miRNA_file, 'w')
        newfile.write('gene' + '\t' + 'transcript_ID' + '\t' + 'Nmirnas' + '\t' + 'Nsites' + '\t' + 'Nsites_per_mirna' + '\t' + 'CNV' + '\n')
        for gene in regulators:
            newfile.write(gene + '\t')
            newfile.write(regulators[gene][-2] + '\t')
            N_mirnas = 0
            N_sites = 0
            for pair in regulators[gene][:-2]:
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/N_mirnas) + '\t' + regulators[gene][-1] + '\n')
        newfile.close()
           
        # get data from file using existing transcript
        # do the test
        # store outcomes in the dictionnary

        # open the temporary outputfile and grab its content
        list_mirnas = get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        
        # get the lists of UTR length
        CNV_UTR, nonCNV_UTR = get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)

        # test that UTR length of CNV and nonCNV genes is different, keep the p-value
        p_UTR = test_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)[1]

        # test that mean number of mirnas, sites and sites per mirna are different netween CNV and nonCNV genes
        test_mirna = test_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        p_mirna, p_sites, p_ratio = test_mirna[0][1], test_mirna[1][1], test_mirna[2][1]

        # compute means for the CNV and non-CNV genes for the 3 variables and for UTR length
        CNV_mean_mirnas = compute_mean_std_error(list_mirnas[0])[0]
        nonCNV_mean_mirnas = compute_mean_std_error(list_mirnas[1])[0]
        CNV_mean_sites = compute_mean_std_error(list_mirnas[2])[0]
        nonCNV_mean_sites = compute_mean_std_error(list_mirnas[3])[0]
        CNV_mean_ratio = compute_mean_std_error(list_mirnas[4])[0]
        nonCNV_mean_ratio = compute_mean_std_error(list_mirnas[5])[0]
        CNV_mean_UTR = compute_mean_std_error(CNV_UTR)[0]
        nonCNV_mean_UTR = compute_mean_std_error(nonCNV_UTR)[0]

        # populate the dictionary
        samples[j] = []
        corrected_samples[j] = []
            
        # compare mean number of mirnas
        if p_mirna < 0.05:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_greater')
        else:
            samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05:
            if CNV_mean_sites < nonCNV_mean_sites:
                samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                samples[j].append('sites_CNV_greater')
        else:
            samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_greater')
        else:
            samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_greater')
        else:
            samples[j].append('UTR_CNV_notdiff')


        # compare mean number of mirnas
        if p_mirna < 0.05 / replicates:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_greater')
        else:
            corrected_samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05 / replicates:
            if CNV_mean_sites < nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_greater')
        else:
            corrected_samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05 / replicates:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_greater')
        else:
            corrected_samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05 / replicates:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_greater')
        else:
            corrected_samples[j].append('UTR_CNV_notdiff')

        # change j and move to the next replicate
        j -= 1
        if j % 5 == 0:
            print(j)
       
    return samples, corrected_samples



def resampling_random_non_CNV_genes(summary_counts, species, conservation_sites, miRNA_family_info, UTR_file, all_gene_file,
                            CNV_file, replicates, ensembl_transcript_coordinates = 'zebrafish_ensembl_genes_zv8',
                            fly_gene_coordinates = 'dmel-all-r5.1.gff', zebrafish_CNV_genes_file = 'Zebrafish_CNV_genes.txt',
                            ensembl_gene_transcripts = 'Zebrafish_Gene_Transcripts_ID', CNV_miRNA_file = 'temp_outputfile.txt'):
    '''
    (file, str, str, file, file, file, file, str, int, file, file, file, file, file) -> dict
    Sample non-CNV genes at random and assign randomly the non-CNV genes to CNV and non-CNV groups. Test for each replicate differences in mean number of mirnas,
    sites and sites per mirna and UTR length between "CNV genes" and non-CNV genes.
    Returns a dictionnary with study as key and a list of outcomes as value (significantly lower, greater, not_diff)
    '''

    import random
    
    # create dictionnary to store the outcomes of the tests for each replicate
    samples = {}
    corrected_samples = {}

    # get the miRNA info for each transcript
    if species == 'worm':
        regulated_transcripts = worm_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)
    elif species == 'human':
        regulated_transcripts = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)
    elif species == 'zebrafish':
        regulated_transcripts = zebrafish_miRNA_regulation_TargetScan_sites(summary_counts, species, miRNA_family_info)
    elif species == 'fly':
        regulated_transcripts = fly_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)

    # get the UTR sequences of all transcripts for each transcript
    UTR_sequences = UTR_sequence(UTR_file, species)
    
    # get the gene ID for each target transcript
    if species == 'fly':
        transcripts = fly_transcript_gene_pairs(UTR_file)
    elif species == 'worm' or species == 'human':
        transcripts = transcript_gene_pairs(summary_counts)
    elif species == 'zebrafish':
        transcripts = {}
        for transcript_ID in regulated_transcripts:
            gene_ID = transcript_ID[:transcript_ID.index('.')] # the transcript ID is not a valid ensembl ID but the gene_ID + .number
            transcripts[transcript_ID] = gene_ID      

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

    # make a reverse dictionnary with gene: [list of transcript]
    target_genes = {}
    for transcript_ID in transcripts:
        gene_ID = transcripts[transcript_ID]
        if gene_ID in target_genes:
            target_genes[gene_ID].append(transcript_ID)
        else:
            target_genes[gene_ID] = [transcript_ID]
        
    # get the set of valid genes
    if species == 'worm':
        valid_genes = sort_valid_worm_genes(all_gene_file)
    elif species == 'human':
        valid_genes = sort_valid_human_genes(all_gene_file)
    elif species == 'zebrafish':
        valid_genes = sort_valid_zebrafish_genes(ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only = 'yes')
    elif species == 'fly':
        valid_genes = sort_valid_fly_genes(fly_gene_coordinates)

    # remove genes that are not in the reference genome
    to_remove = []
    for gene_ID in target_genes:
        if gene_ID not in valid_genes:
            to_remove.append(gene_ID)
    if len(to_remove) != 0:
        for gene_ID in to_remove:
            del target_genes[gene_ID]

    # get the CNV genes
    if species == 'worm':
        CNV_genes = worm_CNV_genes(CNV_file, all_gene_file)
    elif species == 'human':
        CNV_genes = human_CNV_genes(CNV_file, all_gene_file)
    elif species == 'zebrafish':
        CNV_file = zebrafish_CNV_genes_file
        cnv = open(CNV_file, 'r')
        CNV_genes = set()
        for line in cnv:
            line = line.rstrip()
            if line != '':
                CNV_genes.add(line)
        cnv.close()
    elif species == 'fly':
        CNV_genes = set()
        cnv = open(CNV_file, 'r')
        cnv.readline()
        for line in cnv:
            line = line.rstrip()
            if line != '':
                line = line.split()
                CNV_genes.add(line[0])
                CNV_genes.add(line[1])
        cnv.close()
        
    # make a dictionnary of non-CNV genes:
    nonCNV_targets = {}

    for gene in target_genes:
        if gene not in CNV_genes:
            nonCNV_targets[gene] = target_genes[gene]

    # make a list of non-CNV genes
    draw_from_nonCNV_targets = [gene for gene in nonCNV_targets]
    print(len(draw_from_nonCNV_targets))
    

    # repeat the tests for number of replicates, and store the outcomes in the dictionnary of replicates outcomes
    j = replicates
    while j != 0:
        # create a dictionnary to store the miRNA content of each transcript (1 transcript per gene) for CNV genes and for non-CNV genes
        mirna_target_nonCNV = {}
        
        # create a dictionnary to store mirna information for all 500 CNV genes and all 500 non-CNV genes
        regulators = {}
        
        # create a set of genes that are already picked
        already_picked = set()
        
        # pick 1000 non-CNV genes at random and assign CNV at random with 50% probability
        i = 1000
        while i != 0:
            # pick a random number and get the gene from the list of non_cnv_targets using the random index
            position = random.randint(0, len(draw_from_nonCNV_targets)-1)
            target = draw_from_nonCNV_targets[position]
            # if the gene has multiple transcripts, pick a random transcript
            if len(nonCNV_targets[target]) > 1:
                position = random.randint(0, len(nonCNV_targets[target]) -1)
                target_transcript = nonCNV_targets[target][position]
            else:
                target_transcript = nonCNV_targets[target][0]
            
            # verify that the genes were not already picked
            if target not in already_picked:
                # add the genes to set of already picked genes
                already_picked.add(target)
                #and that the transcript is in regulated_transcripts (ie. celegans includes regulation for miR* that were removed)
                if target_transcript in regulated_transcripts:
                    # use the dictionnary of transcript:[mirna info] to extract the mirna info
                    #WARNING: use list to grab the list so that regulated_transcripts is not modified
                    mirna_target_nonCNV[target_transcript] = list(regulated_transcripts[target_transcript]) 
                    mirna_target_nonCNV[target_transcript].append(target_transcript)
                    # pick a random number and use this number to assign CNV status
                    num = random.randint(1, 2)
                    if num % 2 == 0:
                        mirna_target_nonCNV[target_transcript].append('CNV')
                    else:
                        mirna_target_nonCNV[target_transcript].append('non-CNV')
                    # change i value
                    i -= 1
                
        # populate the dictionnary regulators with mirna information, transcript_ID and CNV status
        # write content to temporary file       
        for target_transcript in mirna_target_nonCNV:
            gene = transcripts[target_transcript]
            regulators[gene] = mirna_target_nonCNV[target_transcript]
               
        newfile = open(CNV_miRNA_file, 'w')
        newfile.write('gene' + '\t' + 'transcript_ID' + '\t' + 'Nmirnas' + '\t' + 'Nsites' + '\t' + 'Nsites_per_mirna' + '\t' + 'CNV' + '\n')
        for gene in regulators:
            newfile.write(gene + '\t')
            newfile.write(regulators[gene][-2] + '\t')
            N_mirnas = 0
            N_sites = 0
            for pair in regulators[gene][:-2]:
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/N_mirnas) + '\t' + regulators[gene][-1] + '\n')
        newfile.close()
           
        # get data from file using existing transcript
        # do the test
        # store outcomes in the dictionnary

        # open the temporary outputfile and grab its content
        list_mirnas = get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        
        # get the lists of UTR length
        CNV_UTR, nonCNV_UTR = get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)

        # test that UTR length of CNV and nonCNV genes is different, keep the p-value
        p_UTR = test_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)[1]

        # test that mean number of mirnas, sites and sites per mirna are different netween CNV and nonCNV genes
        test_mirna = test_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        p_mirna, p_sites, p_ratio = test_mirna[0][1], test_mirna[1][1], test_mirna[2][1]

        # compute means for the CNV and non-CNV genes for the 3 variables and for UTR length
        CNV_mean_mirnas = compute_mean_std_error(list_mirnas[0])[0]
        nonCNV_mean_mirnas = compute_mean_std_error(list_mirnas[1])[0]
        CNV_mean_sites = compute_mean_std_error(list_mirnas[2])[0]
        nonCNV_mean_sites = compute_mean_std_error(list_mirnas[3])[0]
        CNV_mean_ratio = compute_mean_std_error(list_mirnas[4])[0]
        nonCNV_mean_ratio = compute_mean_std_error(list_mirnas[5])[0]
        CNV_mean_UTR = compute_mean_std_error(CNV_UTR)[0]
        nonCNV_mean_UTR = compute_mean_std_error(nonCNV_UTR)[0]

        # populate the dictionary
        samples[j] = []
        corrected_samples[j] = []
            
        # compare mean number of mirnas
        if p_mirna < 0.05:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_greater')
        else:
            samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05:
            if CNV_mean_sites < nonCNV_mean_sites:
                samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                samples[j].append('sites_CNV_greater')
        else:
            samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_greater')
        else:
            samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_greater')
        else:
            samples[j].append('UTR_CNV_notdiff')


        # compare mean number of mirnas
        if p_mirna < 0.05 / replicates:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_greater')
        else:
            corrected_samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05 / replicates:
            if CNV_mean_sites < nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_greater')
        else:
            corrected_samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05 / replicates:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_greater')
        else:
            corrected_samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05 / replicates:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_greater')
        else:
            corrected_samples[j].append('UTR_CNV_notdiff')

        # change j and move to the next replicate
        j -= 1
        if j % 5 == 0:
            print(j)
       
    return samples, corrected_samples






def resampling_random_CNV_genes(summary_counts, species, conservation_sites, miRNA_family_info, UTR_file, all_gene_file,
                            CNV_file, replicates, ensembl_transcript_coordinates = 'zebrafish_ensembl_genes_zv8',
                            fly_gene_coordinates = 'dmel-all-r5.1.gff', zebrafish_CNV_genes_file = 'Zebrafish_CNV_genes.txt',
                            ensembl_gene_transcripts = 'Zebrafish_Gene_Transcripts_ID', CNV_miRNA_file = 'temp_outputfile.txt'):
    '''
    (file, str, str, file, file, file, file, str, int, file, file, file, file, file) -> dict
    Sample CNV genes at random and assign randomly the non-CNV genes to CNV and non-CNV groups. Test for each replicate differences in mean number of mirnas,
    sites and sites per mirna and UTR length between "CNV genes" and non-CNV genes.
    Returns a dictionnary with study as key and a list of outcomes as value (significantly lower, greater, not_diff)
    '''

    import random
    
    # create dictionnary to store the outcomes of the tests for each replicate
    samples = {}
    corrected_samples = {}

    # get the miRNA info for each transcript
    if species == 'worm':
        regulated_transcripts = worm_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)
    elif species == 'human':
        regulated_transcripts = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)
    elif species == 'zebrafish':
        regulated_transcripts = zebrafish_miRNA_regulation_TargetScan_sites(summary_counts, species, miRNA_family_info)
    elif species == 'fly':
        regulated_transcripts = fly_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)

    # get the UTR sequences of all transcripts for each transcript
    UTR_sequences = UTR_sequence(UTR_file, species)
    
    # get the gene ID for each target transcript
    if species == 'fly':
        transcripts = fly_transcript_gene_pairs(UTR_file)
    elif species == 'worm' or species == 'human':
        transcripts = transcript_gene_pairs(summary_counts)
    elif species == 'zebrafish':
        transcripts = {}
        for transcript_ID in regulated_transcripts:
            gene_ID = transcript_ID[:transcript_ID.index('.')] # the transcript ID is not a valid ensembl ID but the gene_ID + .number
            transcripts[transcript_ID] = gene_ID      

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

    # make a reverse dictionnary with gene: [list of transcript]
    target_genes = {}
    for transcript_ID in transcripts:
        gene_ID = transcripts[transcript_ID]
        if gene_ID in target_genes:
            target_genes[gene_ID].append(transcript_ID)
        else:
            target_genes[gene_ID] = [transcript_ID]
        
    # get the set of valid genes
    if species == 'worm':
        valid_genes = sort_valid_worm_genes(all_gene_file)
    elif species == 'human':
        valid_genes = sort_valid_human_genes(all_gene_file)
    elif species == 'zebrafish':
        valid_genes = sort_valid_zebrafish_genes(ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only = 'yes')
    elif species == 'fly':
        valid_genes = sort_valid_fly_genes(fly_gene_coordinates)

    # remove genes that are not in the reference genome
    to_remove = []
    for gene_ID in target_genes:
        if gene_ID not in valid_genes:
            to_remove.append(gene_ID)
    if len(to_remove) != 0:
        for gene_ID in to_remove:
            del target_genes[gene_ID]

    # get the CNV genes
    if species == 'worm':
        CNV_genes = worm_CNV_genes(CNV_file, all_gene_file)
    elif species == 'human':
        CNV_genes = human_CNV_genes(CNV_file, all_gene_file)
    elif species == 'zebrafish':
        CNV_file = zebrafish_CNV_genes_file
        cnv = open(CNV_file, 'r')
        CNV_genes = set()
        for line in cnv:
            line = line.rstrip()
            if line != '':
                CNV_genes.add(line)
        cnv.close()
    elif species == 'fly':
        CNV_genes = set()
        cnv = open(CNV_file, 'r')
        cnv.readline()
        for line in cnv:
            line = line.rstrip()
            if line != '':
                line = line.split()
                CNV_genes.add(line[0])
                CNV_genes.add(line[1])
        cnv.close()
        
    # make a dictionnary of non-CNV genes:
    CNV_targets = {}

    for gene in target_genes:
        if gene in CNV_genes:
            CNV_targets[gene] = target_genes[gene]

    # make a list of non-CNV genes
    draw_from_CNV_targets = [gene for gene in CNV_targets]
    print(len(draw_from_CNV_targets))
    

    # repeat the tests for number of replicates, and store the outcomes in the dictionnary of replicates outcomes
    j = replicates
    while j != 0:
        # create a dictionnary to store the miRNA content of each transcript (1 transcript per gene) for CNV genes and for non-CNV genes
        mirna_target_CNV = {}
        
        # create a dictionnary to store mirna information for all 500 CNV genes and all 500 non-CNV genes
        regulators = {}
        
        # create a set of genes that are already picked
        already_picked = set()
        
        # pick 1000 non-CNV genes at random and assign CNV at random with 50% probability
        i = 1000
        while i != 0:
            # pick a random number and get the gene from the list of non_cnv_targets using the random index
            position = random.randint(0, len(draw_from_CNV_targets)-1)
            target = draw_from_CNV_targets[position]
            # if the gene has multiple transcripts, pick a random transcript
            if len(CNV_targets[target]) > 1:
                position = random.randint(0, len(CNV_targets[target]) -1)
                target_transcript = CNV_targets[target][position]
            else:
                target_transcript = CNV_targets[target][0]
            
            # verify that the genes were not already picked
            if target not in already_picked:
                # add the genes to set of already picked genes
                already_picked.add(target)
                #and that the transcript is in regulated_transcripts (ie. celegans includes regulation for miR* that were removed)
                if target_transcript in regulated_transcripts:
                    # use the dictionnary of transcript:[mirna info] to extract the mirna info
                    #WARNING: use list to grab the list so that regulated_transcripts is not modified
                    mirna_target_CNV[target_transcript] = list(regulated_transcripts[target_transcript]) 
                    mirna_target_CNV[target_transcript].append(target_transcript)
                    # pick a random number and use this number to assign CNV status
                    num = random.randint(1, 2)
                    if num % 2 == 0:
                        mirna_target_CNV[target_transcript].append('CNV')
                    else:
                        mirna_target_CNV[target_transcript].append('non-CNV')
                    # change i value
                    i -= 1
                
        # populate the dictionnary regulators with mirna information, transcript_ID and CNV status
        # write content to temporary file       
        for target_transcript in mirna_target_CNV:
            gene = transcripts[target_transcript]
            regulators[gene] = mirna_target_CNV[target_transcript]
               
        newfile = open(CNV_miRNA_file, 'w')
        newfile.write('gene' + '\t' + 'transcript_ID' + '\t' + 'Nmirnas' + '\t' + 'Nsites' + '\t' + 'Nsites_per_mirna' + '\t' + 'CNV' + '\n')
        for gene in regulators:
            newfile.write(gene + '\t')
            newfile.write(regulators[gene][-2] + '\t')
            N_mirnas = 0
            N_sites = 0
            for pair in regulators[gene][:-2]:
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/N_mirnas) + '\t' + regulators[gene][-1] + '\n')
        newfile.close()
           
        # get data from file using existing transcript
        # do the test
        # store outcomes in the dictionnary

        # open the temporary outputfile and grab its content
        list_mirnas = get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        
        # get the lists of UTR length
        CNV_UTR, nonCNV_UTR = get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)

        # test that UTR length of CNV and nonCNV genes is different, keep the p-value
        p_UTR = test_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)[1]

        # test that mean number of mirnas, sites and sites per mirna are different netween CNV and nonCNV genes
        test_mirna = test_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        p_mirna, p_sites, p_ratio = test_mirna[0][1], test_mirna[1][1], test_mirna[2][1]

        # compute means for the CNV and non-CNV genes for the 3 variables and for UTR length
        CNV_mean_mirnas = compute_mean_std_error(list_mirnas[0])[0]
        nonCNV_mean_mirnas = compute_mean_std_error(list_mirnas[1])[0]
        CNV_mean_sites = compute_mean_std_error(list_mirnas[2])[0]
        nonCNV_mean_sites = compute_mean_std_error(list_mirnas[3])[0]
        CNV_mean_ratio = compute_mean_std_error(list_mirnas[4])[0]
        nonCNV_mean_ratio = compute_mean_std_error(list_mirnas[5])[0]
        CNV_mean_UTR = compute_mean_std_error(CNV_UTR)[0]
        nonCNV_mean_UTR = compute_mean_std_error(nonCNV_UTR)[0]

        # populate the dictionary
        samples[j] = []
        corrected_samples[j] = []
            
        # compare mean number of mirnas
        if p_mirna < 0.05:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_greater')
        else:
            samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05:
            if CNV_mean_sites < nonCNV_mean_sites:
                samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                samples[j].append('sites_CNV_greater')
        else:
            samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_greater')
        else:
            samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_greater')
        else:
            samples[j].append('UTR_CNV_notdiff')


        # compare mean number of mirnas
        if p_mirna < 0.05 / replicates:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_greater')
        else:
            corrected_samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05 / replicates:
            if CNV_mean_sites < nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_greater')
        else:
            corrected_samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05 / replicates:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_greater')
        else:
            corrected_samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05 / replicates:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_greater')
        else:
            corrected_samples[j].append('UTR_CNV_notdiff')

        # change j and move to the next replicate
        j -= 1
        if j % 5 == 0:
            print(j)
       
    return samples, corrected_samples




def assigning_random_CNV_status(summary_counts, species, conservation_sites, miRNA_family_info, UTR_file, all_gene_file,
                            CNV_file, replicates, ensembl_transcript_coordinates = 'zebrafish_ensembl_genes_zv8',
                            fly_gene_coordinates = 'dmel-all-r5.1.gff', zebrafish_CNV_genes_file = 'Zebrafish_CNV_genes.txt',
                            ensembl_gene_transcripts = 'Zebrafish_Gene_Transcripts_ID', CNV_miRNA_file = 'temp_outputfile.txt'):
    '''
    (file, str, str, file, file, file, file, str, int, file, file, file, file, file) -> dict
    Sample genes at random and assign randomly the genes to CNV and non-CNV groups. Test for each replicate differences in mean number of mirnas,
    sites and sites per mirna and UTR length between "CNV genes" and non-CNV genes.
    Returns a dictionnary with study as key and a list of outcomes as value (significantly lower, greater, not_diff)
    '''

    import random
    
    # create dictionnary to store the outcomes of the tests for each replicate
    samples = {}
    corrected_samples = {}

    # get the miRNA info for each transcript
    if species == 'worm':
        regulated_transcripts = worm_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)
    elif species == 'human':
        regulated_transcripts = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)
    elif species == 'zebrafish':
        regulated_transcripts = zebrafish_miRNA_regulation_TargetScan_sites(summary_counts, species, miRNA_family_info)
    elif species == 'fly':
        regulated_transcripts = fly_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)

    # get the UTR sequences of all transcripts for each transcript
    UTR_sequences = UTR_sequence(UTR_file, species)
    
    # get the gene ID for each target transcript
    if species == 'fly':
        transcripts = fly_transcript_gene_pairs(UTR_file)
    elif species == 'worm' or species == 'human':
        transcripts = transcript_gene_pairs(summary_counts)
    elif species == 'zebrafish':
        transcripts = {}
        for transcript_ID in regulated_transcripts:
            gene_ID = transcript_ID[:transcript_ID.index('.')] # the transcript ID is not a valid ensembl ID but the gene_ID + .number
            transcripts[transcript_ID] = gene_ID      

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

    # make a reverse dictionnary with gene: [list of transcript]
    target_genes = {}
    for transcript_ID in transcripts:
        gene_ID = transcripts[transcript_ID]
        if gene_ID in target_genes:
            target_genes[gene_ID].append(transcript_ID)
        else:
            target_genes[gene_ID] = [transcript_ID]
        
    # get the set of valid genes
    if species == 'worm':
        valid_genes = sort_valid_worm_genes(all_gene_file)
    elif species == 'human':
        valid_genes = sort_valid_human_genes(all_gene_file)
    elif species == 'zebrafish':
        valid_genes = sort_valid_zebrafish_genes(ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only = 'yes')
    elif species == 'fly':
        valid_genes = sort_valid_fly_genes(fly_gene_coordinates)

    # remove genes that are not in the reference genome
    to_remove = []
    for gene_ID in target_genes:
        if gene_ID not in valid_genes:
            to_remove.append(gene_ID)
    if len(to_remove) != 0:
        for gene_ID in to_remove:
            del target_genes[gene_ID]

    # make a list of genes
    draw_from_targets = [gene for gene in target_genes]
    print(len(draw_from_targets))
    

    # repeat the tests for number of replicates, and store the outcomes in the dictionnary of replicates outcomes
    j = replicates
    while j != 0:
        # create a dictionnary to store the miRNA content of each transcript (1 transcript per gene) for CNV genes and for non-CNV genes
        mirna_target_CNV = {}
        
        # create a dictionnary to store mirna information for all 500 CNV genes and all 500 non-CNV genes
        regulators = {}
        
        # create a set of genes that are already picked
        already_picked = set()
        
        # pick 1000 genes at random and assign CNV at random with 50% probability
        i = 1000
        while i != 0:
            # pick a random number and get the gene from the list of non_cnv_targets using the random index
            position = random.randint(0, len(draw_from_targets)-1)
            target = draw_from_targets[position]
            # if the gene has multiple transcripts, pick a random transcript
            if len(target_genes[target]) > 1:
                position = random.randint(0, len(target_genes[target]) -1)
                target_transcript = target_genes[target][position]
            else:
                target_transcript = target_genes[target][0]
            
            # verify that the genes were not already picked
            if target not in already_picked:
                # add the genes to set of already picked genes
                already_picked.add(target)
                #and that the transcript is in regulated_transcripts (ie. celegans includes regulation for miR* that were removed)
                if target_transcript in regulated_transcripts:
                    # use the dictionnary of transcript:[mirna info] to extract the mirna info
                    #WARNING: use list to grab the list so that regulated_transcripts is not modified
                    mirna_target_CNV[target_transcript] = list(regulated_transcripts[target_transcript]) 
                    mirna_target_CNV[target_transcript].append(target_transcript)
                    # pick a random number and use this number to assign CNV status
                    num = random.randint(1, 2)
                    if num % 2 == 0:
                        mirna_target_CNV[target_transcript].append('CNV')
                    else:
                        mirna_target_CNV[target_transcript].append('non-CNV')
                    # change i value
                    i -= 1
                
        # populate the dictionnary regulators with mirna information, transcript_ID and CNV status
        # write content to temporary file       
        for target_transcript in mirna_target_CNV:
            gene = transcripts[target_transcript]
            regulators[gene] = mirna_target_CNV[target_transcript]
               
        newfile = open(CNV_miRNA_file, 'w')
        newfile.write('gene' + '\t' + 'transcript_ID' + '\t' + 'Nmirnas' + '\t' + 'Nsites' + '\t' + 'Nsites_per_mirna' + '\t' + 'CNV' + '\n')
        for gene in regulators:
            newfile.write(gene + '\t')
            newfile.write(regulators[gene][-2] + '\t')
            N_mirnas = 0
            N_sites = 0
            for pair in regulators[gene][:-2]:
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/N_mirnas) + '\t' + regulators[gene][-1] + '\n')
        newfile.close()
           
        # get data from file using existing transcript
        # do the test
        # store outcomes in the dictionnary

        # open the temporary outputfile and grab its content
        list_mirnas = get_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        
        # get the lists of UTR length
        CNV_UTR, nonCNV_UTR = get_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)

        # test that UTR length of CNV and nonCNV genes is different, keep the p-value
        p_UTR = test_UTR_length_CNV_non_CNV_genes(UTR_file, species, CNV_miRNA_file)[1]

        # test that mean number of mirnas, sites and sites per mirna are different netween CNV and nonCNV genes
        test_mirna = test_miRNA_regulation_CNV_non_CNV_genes(CNV_miRNA_file)
        p_mirna, p_sites, p_ratio = test_mirna[0][1], test_mirna[1][1], test_mirna[2][1]

        # compute means for the CNV and non-CNV genes for the 3 variables and for UTR length
        CNV_mean_mirnas = compute_mean_std_error(list_mirnas[0])[0]
        nonCNV_mean_mirnas = compute_mean_std_error(list_mirnas[1])[0]
        CNV_mean_sites = compute_mean_std_error(list_mirnas[2])[0]
        nonCNV_mean_sites = compute_mean_std_error(list_mirnas[3])[0]
        CNV_mean_ratio = compute_mean_std_error(list_mirnas[4])[0]
        nonCNV_mean_ratio = compute_mean_std_error(list_mirnas[5])[0]
        CNV_mean_UTR = compute_mean_std_error(CNV_UTR)[0]
        nonCNV_mean_UTR = compute_mean_std_error(nonCNV_UTR)[0]

        # populate the dictionary
        samples[j] = []
        corrected_samples[j] = []
            
        # compare mean number of mirnas
        if p_mirna < 0.05:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                samples[j].append('mirnas_CNV_greater')
        else:
            samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05:
            if CNV_mean_sites < nonCNV_mean_sites:
                samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                samples[j].append('sites_CNV_greater')
        else:
            samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                samples[j].append('ratio_CNV_greater')
        else:
            samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                samples[j].append('UTR_CNV_greater')
        else:
            samples[j].append('UTR_CNV_notdiff')


        # compare mean number of mirnas
        if p_mirna < 0.05 / replicates:
            if CNV_mean_mirnas < nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_lower')
            elif CNV_mean_mirnas > nonCNV_mean_mirnas:
                corrected_samples[j].append('mirnas_CNV_greater')
        else:
            corrected_samples[j].append('mirnas_CNV_notdiff')
        # compare mean number of sites
        if p_sites < 0.05 / replicates:
            if CNV_mean_sites < nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_lower')
            elif CNV_mean_sites > nonCNV_mean_sites:
                corrected_samples[j].append('sites_CNV_greater')
        else:
            corrected_samples[j].append('sites_CNV_notdiff')
        # compare mean number of sites per mirna
        if p_ratio < 0.05 / replicates:
            if CNV_mean_ratio < nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_lower')
            elif CNV_mean_ratio > nonCNV_mean_ratio:
                corrected_samples[j].append('ratio_CNV_greater')
        else:
            corrected_samples[j].append('ratio_CNV_notdiff')

        # compare mean UTR length
        if p_UTR < 0.05 / replicates:
            if CNV_mean_UTR < nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_lower')
            elif CNV_mean_UTR > nonCNV_mean_UTR:
                corrected_samples[j].append('UTR_CNV_greater')
        else:
            corrected_samples[j].append('UTR_CNV_notdiff')

        # change j and move to the next replicate
        j -= 1
        if j % 5 == 0:
            print(j)
       
    return samples, corrected_samples











def cel_gene_names(gene_name_file):
    '''
    (file) -> dict
    Returns a dictionnary with the WBGene name as key and the C. elegans gene names and gene IDs as value
    (eg: WBGene00000278:[cad-1], ,WBGene00000007: [aat-6,T11F9.4])
    '''
    gene_names = {}
    cel = open(gene_name_file, 'r')
    for line in cel:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            while '' in line:
                line.remove('')
            if line[-1] == 'Live':
                if len(line) == 4:
                    gene_names[line[1]] = [line[1]]
                elif len(line) == 5:
                    gene_names[line[1]] = [line[2], line[3]]
                
    cel.close()
    return gene_names


def cel_CNV_type(CNV_file, gene_name_file):
    '''
    (file, file) -> tuple
    Returns a tuple with 2 dictionnaries corresponding to genes affected by amplification CNV and to genes affected by deletion CNV
    Each dictionnary is in the form {WBGene1: [name1, name2]}, WBgene2:[name]}
    '''
    # create a dictionnary of CNV type
    CNV_type = {}

    # go through the CNV file and get the type of CNV for all the CNV genes
    CNV = open(CNV_file, 'r')
    CNV.readline()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            line = line.split()
            CNV_type[line[2]] = line[6]

    #get the gene names
    gene_names = cel_gene_names(gene_name_file)

    # make a set of CNV-deletion and CNV-amplification genes
    amplification = {}
    deletion = {}
    for gene in CNV_type:
        if CNV_type[gene] == 'A':
            amplification[gene] = CNV_type[gene]
        elif CNV_type[gene] == 'D':
            deletion[gene] = CNV_type[gene]

    CNV.close()
    return amplification, deletion



def cel_get_miRNA_regulation_CNV_type(CNV_miRNA_file, CNV_file, gene_name_file):
    '''
    (file) -> tuple
    Returns a 9-item tuple containing the lists of the number of miRNA regulators,
    the number of miRNA binding sites, and the number of sites per miRN for CNV-amplification, CNV-deletion and non-CNV genes
    '''

    # get the set of genes for amplification and deletion CNVs
    amplification, deletion = cel_CNV_type(CNV_file, gene_name_file)
    
    # open CNV_miRNA file
    CNV_miRNA = open(CNV_miRNA_file, 'r')
    
    # create lists to store the values of the different metrics for CNV-Amplification, CNV-deletion and non-CNV genes
    CNV_Del_N_mirnas = []
    CNV_Del_N_sites = []
    CNV_Del_ratio = []
    CNV_Ampl_N_mirnas = []
    CNV_Ampl_N_sites = []
    CNV_Ampl_ratio = []
    nonCNV_N_mirnas = []
    nonCNV_N_sites = []
    nonCNV_ratio = []

    # read file and store the the values in the lists
    CNV_miRNA.readline()
    for line in CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[-1] == 'CNV':
                for gene in amplification:
                    if line[0] in amplification[gene]:
                        CNV_Ampl_N_mirnas.append(int(line[2]))
                        CNV_Ampl_N_sites.append(int(line[3]))
                        CNV_Ampl_ratio.append(float(line[4]))
                for gene in deletion:
                    if line[0] in deletion[gene]:
                        CNV_Del_N_mirnas.append(int(line[2]))
                        CNV_Del_N_sites.append(int(line[3]))
                        CNV_Del_ratio.append(float(line[4]))
            elif line[-1] == 'non-CNV':
                nonCNV_N_mirnas.append(int(line[2]))
                nonCNV_N_sites.append(int(line[3]))
                nonCNV_ratio.append(float(line[4]))
        
    CNV_miRNA.close()
    return CNV_Ampl_N_mirnas, CNV_Del_N_mirnas, nonCNV_N_mirnas, CNV_Ampl_N_sites, CNV_Del_N_sites, nonCNV_N_sites, CNV_Ampl_ratio, CNV_Del_ratio, nonCNV_ratio




def cel_test_miRNA_regulation_CNV_type(CNV_miRNA_file, CNV_file, gene_name_file):
    '''
    (file) -> tuple
    Performs a Wilcoxon rank sum test between CNV-amplification, CNV-deletion and non-CNV genes for
    the number of miRNA regulators, the number of miRNA binding sites, and the number of sites per miRNA
    Return a tuple containing 6 2-item tuples with the z-value and the p-value
    '''

    # get the number of miRNAs, sites and sites per miRNA for CNV and non-CNV genes
    miRNA_regulation = cel_get_miRNA_regulation_CNV_type(CNV_miRNA_file, CNV_file, gene_name_file)
    CNV_Ampl_N_mirnas = miRNA_regulation[0]
    CNV_Del_N_mirnas = miRNA_regulation[1]
    nonCNV_N_mirnas = miRNA_regulation[2]
    CNV_Ampl_N_sites = miRNA_regulation[3]
    CNV_Del_N_sites = miRNA_regulation[4]
    nonCNV_N_sites = miRNA_regulation[5]
    CNV_Ampl_ratio = miRNA_regulation[6]
    CNV_Del_ratio = miRNA_regulation[7]
    nonCNV_ratio = miRNA_regulation[8]

    # compute the Wilcoxon rank sum tests
    from scipy import stats
    diff_N_mirnas_ampl_del = stats.ranksums(CNV_Ampl_N_mirnas, CNV_Del_N_mirnas)
    diff_N_sites_ampl_del = stats.ranksums(CNV_Ampl_N_sites, CNV_Del_N_sites)
    diff_ratio_ampl_del = stats.ranksums(CNV_Ampl_ratio, CNV_Del_ratio)
    diff_N_mirnas_ampl_noncnv = stats.ranksums(CNV_Ampl_N_mirnas, nonCNV_N_mirnas)
    diff_N_sites_ampl_noncnv = stats.ranksums(CNV_Ampl_N_sites, nonCNV_N_sites)
    diff_ratio_ampl_noncnv = stats.ranksums(CNV_Ampl_ratio, nonCNV_ratio)
    diff_N_mirnas_del_noncnv = stats.ranksums(CNV_Del_N_mirnas, nonCNV_N_mirnas)
    diff_N_sites_del_noncnv = stats.ranksums(CNV_Del_N_sites, nonCNV_N_sites)
    diff_ratio_del_noncnv = stats.ranksums(CNV_Del_ratio, nonCNV_ratio)

    return ((diff_N_mirnas_ampl_del, diff_N_mirnas_ampl_noncnv, diff_N_mirnas_del_noncnv),
            (diff_N_sites_ampl_del, diff_N_sites_ampl_noncnv, diff_N_sites_del_noncnv), 
           (diff_ratio_ampl_del, diff_ratio_ampl_noncnv, diff_ratio_del_noncnv))
