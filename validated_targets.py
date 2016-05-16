from CNV_miRNAs_TargetScan import *
from CNV_nonCNV_genes_comparisons import *
from scipy import stats


def get_mirna_interactions_miRecords(validated_targets_file, species):
    '''
    (file, str) -> dict
    Return a dictionnary with genes as keys and a set of miRNA regulators as value
    '''

    # open target_file for reading
    infile = open(validated_targets_file, 'r')
    header = infile.readline().rstrip()

    # create a dictionary to store the miRNA regulators
    targets = {}

    # go through the file and record the mirnas for each gene
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            gene = line[2]
            mirna = line[-1]
            if species == 'human':
                # keep interactions for genes and mirnas from the same species
                if line[1] == 'Homo_sapiens' and line[-2] == 'Homo_sapiens':
                    if gene in targets:
                        targets[gene].add(line[-1])
                    else:
                        targets[gene] = {line[-1],}
            elif species == 'worm':
                if line[1] == 'Caenorhabditis_elegans' and line[-2] == 'Caenorhabditis_elegans':
                    if gene in targets:
                        targets[gene].add(line[-1])
                    else:
                        targets[gene] = {line[-1],}
            elif species == 'fly':
                if line[1] == 'Drosophila_melanogaster' and line[-2] == 'Drosophila_melanogaster':
                    if gene in targets:
                        targets[gene].add(line[-1])
                    else:
                        targets[gene] = {line[-1],}
            elif species == 'fish':
                if line[1] == 'Danio_rerio' and line[-2] == 'Danio_rerio':
                    if gene in targets:
                        targets[gene].add(line[-1])
                    else:
                        targets[gene] = {line[-1],}
    #close file after reading
    infile.close()

    return targets
        
    

    
def CNV_validated_targets_miRecords(validated_targets_file, CNV_file, all_gene_file, species):
    '''
    (file, file, file, str) -> dict
    Returns a dictionnary of target_genes : [{mirna regulators}, CNV_status] for validated miRNA interactions
    in human. Genes not found in the current genome assembly are removed
    '''

    # get the validated interactions from miRecords
    interactions = get_mirna_interactions_miRecords(validated_targets_file, species)

    # get the set of valid genes
    valid_genes = sort_valid_human_genes(all_gene_file)

    # remove non-valid genes from the validated targets
    non_valid = []
    for gene in interactions:
        if gene not in valid_genes:
            non_valid.append(gene)
    for gene in non_valid:
        del interactions[gene]
    
    # get the set of human CNV genes
    CNV_genes = human_CNV_genes(CNV_file, all_gene_file)

    # create a new dictionnary with CNV status
    CNV_interactions = {}

    # add CNV status to all genes in validated targets
    for gene in interactions:
        CNV_interactions[gene] = []
        CNV_interactions[gene].append(interactions[gene])
        if gene in CNV_genes:
            CNV_interactions[gene].append('CNV')
        else:
            CNV_interactions[gene].append('non-CNV')

    return CNV_interactions


def compare_miRNAs_validated_targets(validated_targets_file, CNV_file, all_gene_file, database, strength, evidence, species):
    '''
    (file, file, file, str, str) -> tuple of tuples
    Performs a Wilcoxon ranksum test for mean difference in the number of
    miRNA regulators between validated CNV and non_CNV target genes
    and return the mean # miRNAs and SEM for CNVs, non-CNV genes and with the z-score and the p-value
    The database source for the validated targets needs to be specified and the strength of interaction can be indicated
    when the database is miRTarbase
    '''

    # get the validated targets with CNV status
    if database == 'miRecords':
        interactions = CNV_validated_targets_miRecords(validated_targets_file, CNV_file, all_gene_file, species)
    elif database == 'miRTarbase':
        interactions = CNV_validated_targets_miRTarbase(validated_targets_file, CNV_file, all_gene_file, strength, evidence, species)
    
    # make lists for CNV and non_CNV genes
    CNV_mirnas = []
    non_CNV_mirnas = []

    # populate the lists
    for gene in interactions:
        if interactions[gene][-1] == 'CNV':
            CNV_mirnas.append(len(interactions[gene][0]))
        elif interactions[gene][-1] == 'non-CNV':
            non_CNV_mirnas.append(len(interactions[gene][0]))

    # compute mean and standard error of the mean (SEM) for the CNV genes
    mean_mirnas_CNV = compute_mean_std_error(CNV_mirnas)

    # compute mean and standard error of the mean (SEM) for the non_CNV genes
    mean_mirnas_non_CNV = compute_mean_std_error(non_CNV_mirnas)

    # perform a Wilcoxon sum rank tests of mean differences
    diff_mirnas = stats.ranksums(CNV_mirnas, non_CNV_mirnas)

    return len(CNV_mirnas), len(non_CNV_mirnas), mean_mirnas_CNV, mean_mirnas_non_CNV, diff_mirnas

    
    
    

    
def get_mirna_interactions_miRTarbase(validated_targets_file, strength, evidence, species):
    '''
    (file, str, str, str) -> dict
    Return a dictionnary with genes as keys and a set of miRNA regulators as value
    for human or worm keeping all interactions of only the strongest and filtering by the type of evididence
    '''

    # open target_file for reading
    infile = open(validated_targets_file, 'r')
    header = infile.readline().rstrip()

    # create a dictionary to store the miRNA regulators
    targets = {}

    # make a set of keywords corresponding to strong evidence (western-blot, reporter assays, qPCR and sequencing)
    strong_evidence = {'Western_blot', 'Reporter_assay', 'reporter_assay',
                       'Luciferase_assay', 'LacZ_assay', 'Real_Time_RT-PCR',
                       'qRT_PCR', 'Real_time_PCR', 'CLASH', 'ChIP', 'Sequencing', 'CLASH', 'CLIP-seq'}
                       

    # go through the file and record the mirnas for each gene
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            gene = line[3]
            mirna = line[1]
            # ignore genes with non-functional interactions
            if 'Non-Functional' not in line[7]:
                # use all the evidence of the strongest
                if evidence == 'all_evidence':
                    # use all interactions or the strongest
                    if strength == 'all':
                        if species == 'human':
                            # keep interactions for genes and mirnas from the same species
                            if line[2] == 'Homo_sapiens' and line[5] == 'Homo_sapiens':
                                if gene in targets:
                                    targets[gene].add(mirna)
                                else:
                                    targets[gene] = {mirna,}
                        elif species == 'worm':
                            # keep interactions for genes and mirnas from the same species
                            if line[2] == 'Caenorhabditis_elegans' and line[5] == 'Caenorhabditis_elegans':
                                if gene in targets:
                                    targets[gene].add(mirna)
                                else:
                                    targets[gene] = {mirna,}
                    elif strength == 'strongest':
                        if 'Weak' not in line[7]:
                            if species == 'human':
                                # keep interactions for genes and mirnas from the same species
                                if line[2] == 'Homo_sapiens' and line[5] == 'Homo_sapiens':
                                    if gene in targets:
                                        targets[gene].add(mirna)
                                    else:
                                        targets[gene] = {mirna,}
                            elif species == 'worm':
                                # keep interactions for genes and mirnas from the same species
                                if line[2] == 'Caenorhabditis_elegans' and line[5] == 'Caenorhabditis_elegans':
                                    if gene in targets:
                                        targets[gene].add(mirna)
                                    else:
                                        targets[gene] = {mirna,}
                elif evidence == 'strong_evidence':
                    # check if the method listed corresponds to a strong evidence
                    for item in strong_evidence:
                        if item in line[6]:
                            # use all interactions or the strongest
                            if strength == 'all':
                                if species == 'human':
                                    # keep interactions for genes and mirnas from the same species
                                    if line[2] == 'Homo_sapiens' and line[5] == 'Homo_sapiens':
                                        if gene in targets:
                                            targets[gene].add(mirna)
                                        else:
                                            targets[gene] = {mirna,}
                                elif species == 'worm':
                                    # keep interactions for genes and mirnas from the same species
                                    if line[2] == 'Caenorhabditis_elegans' and line[5] == 'Caenorhabditis_elegans':
                                        if gene in targets:
                                            targets[gene].add(mirna)
                                        else:
                                            targets[gene] = {mirna,}
                            elif strength == 'strongest':
                                if 'Weak' not in line[7]:
                                    if species == 'human':
                                        # keep interactions for genes and mirnas from the same species
                                        if line[2] == 'Homo_sapiens' and line[5] == 'Homo_sapiens':
                                            if gene in targets:
                                                targets[gene].add(mirna)
                                            else:
                                                targets[gene] = {mirna,}
                                    elif species == 'worm':
                                        # keep interactions for genes and mirnas from the same species
                                        if line[2] == 'Caenorhabditis_elegans' and line[5] == 'Caenorhabditis_elegans':
                                            if gene in targets:
                                                targets[gene].add(mirna)
                                            else:
                                                targets[gene] = {mirna,}
    # close file after reading
    infile.close()

    return targets


def CNV_validated_targets_miRTarbase(validated_targets_file, CNV_file, all_gene_file, strength, evidence, species):
    '''
    (file, file, file, str, str) -> dict
    Returns a dictionnary of target_genes : [{mirna regulators}, CNV_status] for validated miRNA interactions
    in human or worm, keeping all interactions or only the strongest ones. Genes not found in the current genome assembly are removed
    '''

    # get the validated interactions from miRTarbase
    interactions = get_mirna_interactions_miRTarbase(validated_targets_file, strength, evidence, species)

    # get the set of valid genes
    if species == 'human':
        valid_genes = sort_valid_human_genes(all_gene_file)
    elif species == 'worm':
        valid_genes = sort_valid_worm_genes(all_gene_file)

    # remove non-valid genes from the validated targets
    non_valid = []
    for gene in interactions:
        if gene not in valid_genes:
            non_valid.append(gene)
    for gene in non_valid:
        del interactions[gene]
    
    # get the set of CNV genes
    if species == 'human':
        CNV_genes = human_CNV_genes(CNV_file, all_gene_file)
    elif species == 'worm':
        CNV_genes = worm_CNV_genes(CNV_file, all_gene_file)

    # create a new dictionnary with CNV status
    CNV_interactions = {}

    # add CNV status to all genes in validated targets
    for gene in interactions:
        CNV_interactions[gene] = []
        CNV_interactions[gene].append(interactions[gene])
        if gene in CNV_genes:
            CNV_interactions[gene].append('CNV')
        else:
            CNV_interactions[gene].append('non-CNV')

    return CNV_interactions



def get_UTR_length_CNV_validated(UTR_file, CNV_miRNA_file, validated_targets_file, CNV_file, all_gene_file, strength, evidence, database,  species):
    '''
    (file, file, file, file, file, str, str, str, str) -> dict
    Returns a dictionnary with genes as keys and a list containing the number of validated miRNA regulators , the UTR length and the CNV status
    '''

    # get the sequence of all UTRs for the given species
    all_UTR = UTR_sequence(UTR_file, species)

    # make a dictionnary to store the transcript : gene pairs
    transcripts_genes = {}
    CNV_miRNA = open(CNV_miRNA_file, 'r')
    CNV_miRNA.readline()
    for line in CNV_miRNA:
        line =line.rstrip()
        if line != '':
            line = line.split()
            transcripts_genes[line[1]] = line[0]
    CNV_miRNA.close()
    
    if species == 'worm':
        # make a dictionnary to hold the worm gene names pairs {'abf1' : 'C50F2.9')
        gene_names = {}
        # make a dictionnary to hold the worm gene ID : gene name pairs {'C50F2.9' : 'abf1')
        gene_IDs = {}
        # open file for reading
        IDs = open(all_gene_file, 'r')
        for line in IDs:
            line = line.rstrip()
            if line != '':
                line = line.split(',')
                if line[-1] == 'Live':
                    while '' in line:
                        line.remove('')
                    line.remove('6239')
                    line.remove('Live')
                    for item in line:
                        if item.startswith('WBGene'):
                            line.remove(item)
                    if len(line) > 1:
                        gene_names[line[0]] = line[1]
                        gene_IDs[line[1]] = line[0]
        IDs.close()

    # make a dictionnary to store the gene : UTR length
    gene_UTR = {}
    
    if species == 'human':
        for transcript in transcripts_genes:
            gene_UTR[transcripts_genes[transcript]] = len(all_UTR[transcript])
    elif species  =='worm':
        # make a dictionnary with the transcripts : gene names
        transcripts_gene_names = {}
        # make a dictionnary with the transcripts : gene IDs
        transcripts_gene_IDs = {}
        for transcript in transcripts_genes:
            if transcripts_genes[transcript] in gene_names:
                transcripts_gene_IDs[transcript] = gene_names[transcripts_genes[transcript]]
            elif transcripts_genes[transcript] in gene_IDs:
                transcripts_gene_names[transcript] = gene_IDs[transcripts_genes[transcript]]
        # populate the dictionnary with gene names and gene Ids paired with the UTR length
        for transcript in transcripts_genes:
            gene_UTR[transcripts_genes[transcript]] = len(all_UTR[transcript])
        for transcript in transcripts_gene_names:
            gene_UTR[transcripts_gene_names[transcript]] = len(all_UTR[transcript])
        for transcript in transcripts_gene_IDs:
            gene_UTR[transcripts_gene_IDs[transcript]] = len(all_UTR[transcript])
      
    
    # get the CNV status of the validated targets
    if database == 'miRTarbase':
        CNV_interactions = CNV_validated_targets_miRTarbase(validated_targets_file, CNV_file, all_gene_file, strength, evidence, species)
    elif database == 'miRecords':
        # this function is only designed for human, too few genes for other species
        CNV_interactions = CNV_validated_targets_miRecords(validated_targets_file, CNV_file, all_gene_file, species)


    # create a dictionnary to store the number of miRNAs, the UTR length and the CNV status {gene : [N_mirnas, UTR_length, CNV]}
    CNV_targets = {}
    for gene in CNV_interactions:
        if gene in gene_UTR:
            CNV_targets[gene] = [len(CNV_interactions[gene][0]), gene_UTR[gene], CNV_interactions[gene][1]]

    return CNV_targets

        
def test_UTR_length_CNV_validayed(UTR_file, CNV_miRNA_file, validated_targets_file, CNV_file, all_gene_file, strength, evidence, database,  species):
    '''
    (file, file, file, file, file, str, str, str, str) -> tuples of tuples
    Return the mean number of validated miRNA regulators for CNV genes and non-CNV genes, the z-score and the p-value of the Wilcoxon sum rank test of mean differences,
    and the Spearman's rho and its p-value for the correlation between UTR length and the number of validated regulators
    '''
    
    # get the number of miRNAs and UTR length for each validated targets
    CNV_targets = get_UTR_length_CNV_validated(UTR_file, CNV_miRNA_file, validated_targets_file, CNV_file, all_gene_file, strength, evidence, database,  species)

    # make lists to the store the UTR length of CNV and non-CNV target transcripts
    CNV_UTR = []
    nonCNV_UTR = []

    # make a list to store the UTR length and the number of miRNAs
    UTR = []
    miRNAs = []

    CNV_miRNAs = []
    nonCNV_miRNAs = []

    # partition the UTR length of transcripts target based on CNV status
    # populate the lists of UTR lengths and number of miRNAs in the same order to compute correlation
    for gene in CNV_targets:
        if CNV_targets[gene][-1] == 'CNV':
            CNV_UTR.append(CNV_targets[gene][-2])
            CNV_miRNAs.append(CNV_targets[gene][0])
        elif CNV_targets[gene][-1] == 'non-CNV':
            nonCNV_UTR.append(CNV_targets[gene][-2])
            nonCNV_miRNAs.append(CNV_targets[gene][0])
        UTR.append(CNV_targets[gene][-2])
        miRNAs.append(CNV_targets[gene][0])

    # compute mean UTR lengths
    mean_UTR_CNV = compute_mean_std_error(CNV_UTR)
    mean_UTR_nonCNV = compute_mean_std_error(nonCNV_UTR)

    # compute mean mirnas
    mean_mirnas_CNV = compute_mean_std_error(CNV_miRNAs)
    mean_mirnas_nonCNV = compute_mean_std_error(nonCNV_miRNAs)

    # perform a Wilcoxon rank correlation between the UTR length and the number of miRNAs
    correlation = stats.spearmanr(UTR, miRNAs)

    # perform a Wilcoxon rank sum test between the UTR length and the number of miRNAs
    diff_mirnas = stats.ranksums(CNV_miRNAs, nonCNV_miRNAs)

    # performa Wilcoxon rank sum test of mean difference between CNv and non CNV genes
    diff_UTR = stats.ranksums(CNV_UTR, nonCNV_UTR)

    return (len(CNV_miRNAs), len(nonCNV_miRNAs)), mean_mirnas_CNV, mean_mirnas_nonCNV, diff_mirnas, mean_UTR_CNV, mean_UTR_nonCNV, diff_UTR, correlation

    
