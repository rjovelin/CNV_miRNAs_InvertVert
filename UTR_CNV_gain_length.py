from CNV_miRNAs_TargetScan import UTR_sequence
from CNV_miRNAs_TargetScan import sort_valid_human_genes
from CNV_nonCNV_genes_comparisons import compute_mean_std_error


def worm_CNV_type(CNV_gene_file):
    '''
    (file) -> dict
    Returns a dictionnary of WBGene as key and CNV type as value for all CNV genes in CNV_gene_file
    '''

    # make a dictionnary of WBGene : CNV type
    indel = {}
    cnv = open(CNV_gene_file, 'r')
    cnv.readline()
    for line in cnv:
        if line.rstrip() != '':
            line = line.split()
            indel[line[2]] = line[6]
    cnv.close()
    return indel


def worm_transcript_to_wormbase_gene(UTR_file):
    '''
    (file) -> dict
    Return a dictionnary with the transcript ID as key and the wormbase WBGene ID as value from the UTR_file
    '''

    # make a dictionnary of transcript : WBGene ID
    utr = open(UTR_file, 'r')
    header = utr.readline()
    WBGenes = {}
    for line in utr:
        line = line.rstrip()
        if line != '':
            line = line.split()            
            if '6239' in line:
                transcript = line[0]
                gene = line[3]
                WBGenes[transcript] = gene
    utr.close()
    return WBGenes

def worm_CNV_nonCNV_transcripts(CNV_miRNA_file):
    '''
    (file) -> (set, set)
    Returns a tuple with the set of CNV and non-CNV transcripts from the CNV_miRNA_file
    '''

    # make a set of CNV transcripts using the TargetScan mirna file
    CNVs = set()
    nonCNVs = set()
    cnv = open(CNV_miRNA_file, 'r')
    cnv.readline()
    for line in cnv:
        if line.rstrip() != '':
            line = line.split()
            if line[-1] == 'CNV':
                CNVs.add(line[1])
            elif line[-1] == 'non-CNV':
                nonCNVs.add(line[1])
    cnv.close()
    return CNVs, nonCNVs


def worm_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_gene_file):
    '''
    (file, str, file, file) -> tuple
    Returns a tuple with the lists of UTR length of miRNA target genes from the TargetScan_mirna file
    that in CNV loss, in CNV gain or not in CNV
    '''

    # get the UTR sequence of each transcript
    UTR_seq = UTR_sequence(UTR_file, species)

    # make a dictionnary of transcript : WBGene 
    WBGenes = worm_transcript_to_wormbase_gene(UTR_file)

    # make a dictionnary of WBGene : CNV type
    indel = worm_CNV_type(CNV_gene_file)

    # make sets of CNV and non-CNV transcripts
    CNVs, nonCNVs = worm_CNV_nonCNV_transcripts(CNV_miRNA_file)

    # store the length of UTR for gain CNVs and loss CNVs
    gain_CNV_UTR = []
    loss_CNV_UTR = []
    non_CNV_UTR = []
    for transcript in WBGenes:
        UTR_length = len(UTR_seq[transcript])
        gene = WBGenes[transcript]
        if transcript in CNVs:            
            if indel[gene] == 'A':
                gain_CNV_UTR.append(UTR_length)
            elif indel[gene] == 'D':
                loss_CNV_UTR.append(UTR_length)
        elif transcript in nonCNVs:
            non_CNV_UTR.append(UTR_length)

    return gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR


def worm_mirnas_cnv_gain_loss(UTR_file, CNV_miRNA_file, CNV_gene_file):
    '''
    (file, str, file, file) -> dict
    Returns a dictionnary with transcript as key and a list with the number of mirnas,
    number of sites, ratio and CNV type status from the CNV_miRNA_file
    '''

    # make a dictionnary of transcript : WBGene 
    WBGenes = worm_transcript_to_wormbase_gene(UTR_file)

    # make a dictionnary of WBGene : CNV type
    indel = worm_CNV_type(CNV_gene_file)

    # make a dictionnary of mirna regulation with transcript as key
    cnv_regulation = {}
    cnv_mirnas = open(CNV_miRNA_file, 'r')
    cnv_mirnas.readline()
    for line in cnv_mirnas:
        if line.rstrip() != '':
            line = line.split()
            if line[-1] == 'CNV':
                transcript = line[1]
                gene = WBGenes[transcript]
                if indel[gene] == 'A':
                    cnv_regulation[transcript] = [int(line[2]), int(line[3]), float(line[4]), 'A']
                elif indel[gene] == 'D':
                    cnv_regulation[transcript] = [int(line[2]), int(line[3]), float(line[4]), 'D']
            elif line[-1] == 'non-CNV':
                transcript = line[1]
                cnv_regulation[transcript] = [int(line[2]), int(line[3]), float(line[4]), 'non-CNV']
                
    cnv_mirnas.close()
    return cnv_regulation


def worm_mirna_regulation_CNV_loss_gain_lists(UTR_file, CNV_miRNA_file, CNV_gene_file):
    '''
    (file, file, file) -> tuple
    Returns a tuple with the lists of number of mirans, sites and ratio for genes located in CNV gain, CNV loss and non-CNV
    '''

    # get mirna information for each mirna target in CNV gain, CNV loss and not in CNV
    cnv_regulation = worm_mirnas_cnv_gain_loss(UTR_file, CNV_miRNA_file, CNV_gene_file)

    # make lists of variables for gain and loss CNV
    mirnas_gain_CNV = []    
    mirnas_loss_CNV = []
    mirnas_non_CNV = []
    sites_gain_CNV = []
    sites_loss_CNV = []
    sites_non_CNV = []
    ratio_gain_CNV = []
    ratio_loss_CNV = []
    ratio_non_CNV = []
    
    # populate lists
    for transcript in cnv_regulation:
        if cnv_regulation[transcript][-1] == 'A':
            mirnas_gain_CNV.append(cnv_regulation[transcript][0])
            sites_gain_CNV.append(cnv_regulation[transcript][1])
            ratio_gain_CNV.append(cnv_regulation[transcript][2])
        elif cnv_regulation[transcript][-1] == 'D':
            mirnas_loss_CNV.append(cnv_regulation[transcript][0])
            sites_loss_CNV.append(cnv_regulation[transcript][1])
            ratio_loss_CNV.append(cnv_regulation[transcript][2])
        elif cnv_regulation[transcript][-1] == 'non-CNV':
            mirnas_non_CNV.append(cnv_regulation[transcript][0])
            sites_non_CNV.append(cnv_regulation[transcript][1])
            ratio_non_CNV.append(cnv_regulation[transcript][2])

    return (mirnas_gain_CNV, mirnas_loss_CNV, mirnas_non_CNV,
            sites_gain_CNV, sites_loss_CNV, sites_non_CNV,
            ratio_gain_CNV, ratio_loss_CNV, ratio_non_CNV)


def human_CNV_type(CNV_file, all_gene_file):
    '''
    (file) -> (set, set)
    Returns a tuple containing non-overlapping sets of genes that are in CNV loss and in CNV gain
    '''

    # open the CNV file for reading
    human_CNV = open('GRCh37_hg19_variants_2013-05-31.txt', 'r')
    human_CNV.readline()

    # make sets of CNV loss and CNV gain
    CNV_loss = set()
    CNV_gain = set()

    # go through the file and and the genes in the 2 sets
    for line in human_CNV:
        if line.rstrip() != '':
            line = line.split('\t')
            if line[4] == 'CNV':
                if ',' in line[-2]:
                    genes = line[-2].split(',')
                    if line[5] == 'Loss' or line[5] == 'Deletion':
                        for gene in genes:
                            CNV_loss.add(gene)
                    elif line[5] == 'Gain' or line[5] == 'Insertion' or line[5] == 'Duplication':
                        for gene in genes:
                            CNV_gain.add(gene)
                elif ',' not in line[-2]:
                    if line[5] == 'Loss' or line[5] == 'Deletion':
                        CNV_loss.add(line[-2])
                    elif line[5] == 'Gain' or line[5] == 'Insertion' or line[5] == 'Duplication':
                        CNV_gain.add(line[-2])
                if ',' in line[-1]:
                    genes = line[-1].split(',')
                    if line[5] == 'Loss' or line[5] == 'Deletion':
                        for gene in genes:
                            CNV_loss.add(gene)
                    elif line[5] == 'Gain' or line[5] == 'Insertion' or line[5] == 'Duplication':
                        for gene in genes:
                            CNV_gain.add(gene)
                elif ',' not in line[-1]:
                    if line[5] == 'Loss' or line[5] == 'Deletion':
                        CNV_loss.add(line[-1])
                    elif line[5] == 'Gain' or line[5] == 'Insertion' or line[5] == 'Duplication':                    
                        CNV_gain.add(line[-1])

    # get the set of valid genes
    human_valid_genes = sort_valid_human_genes(all_gene_file)

    # remove non-valid genes from the CNV gene sets
    to_remove = []
    for gene in CNV_gain:
        if gene not in human_valid_genes:
            to_remove.append(gene)
    for gene in CNV_loss:
        if gene not in human_valid_genes:
            to_remove.append(gene)
    for gene in to_remove:
        if gene in CNV_gain:
            CNV_gain.discard(gene)
        if gene in CNV_loss:
            CNV_loss.discard(gene)

    # remove genes that are both in CNV gain and CNV loss
    CNV_complex = CNV_gain.intersection(CNV_loss)
    for gene in CNV_complex:
        CNV_gain.discard(gene)
        CNV_loss.discard(gene)
    
    human_CNV.close()
    return CNV_gain, CNV_loss



def human_transcript_to_gene(CNV_miRNA_file):
    '''
    (file) -> dict
    Returns a dictionnary of transcript : gene pairs from the CNV_miRNA_file
    '''
    # make a dictionnary of the transcript : gene pairs from the CNV_miRNA_file
    transcript_gene = {}
    mirnas = open(CNV_miRNA_file, 'r')
    mirnas.readline()
    for line in mirnas:
        if line.rstrip() != '':
            line = line.split()
            transcript_gene[line[1]] = line[0]
    mirnas.close()
    return transcript_gene


def human_targets_CNV_type(CNV_miRNA_file, CNV_gene_file, all_gene_file):
    '''
    (file) -> (set, set, set)
    Return a 3-item tuple containing the sets of miRNA target transcripts that are exclusively in CNV gain, CNV loss and non-CNV from the file CNV_miRNA_file
    '''

    # make sets of miRNA transcript targets in CNV gain, CNV loss and non-CNV from CNV_miRNA_file
    transcript_CNV_gain = set()
    transcript_CNV_loss = set()
    transcript_non_CNV = set()

    # make a dictionnary of the transcript : gene pairs from the CNV_miRNA_file
    transcript_gene = human_transcript_to_gene(CNV_miRNA_file)

    # get all genes in CNV gain and CNV loss
    CNV_gain, CNV_loss = human_CNV_type(CNV_gene_file, all_gene_file)

    mirnas = open(CNV_miRNA_file, 'r')
    mirnas.readline()
    for line in mirnas:
        if line.rstrip() != '':
            line = line.split()
            transcript = line[1]
            if line[-1] == 'non-CNV':
                transcript_non_CNV.add(transcript)
            elif line[-1] == 'CNV':
                if transcript_gene[transcript] in CNV_gain:
                    transcript_CNV_gain.add(transcript)
                elif transcript_gene[transcript] in CNV_loss:
                    transcript_CNV_loss.add(transcript)
    mirnas.close()
    return transcript_CNV_gain, transcript_CNV_loss, transcript_non_CNV

def human_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_gene_file, all_gene_file):
    '''
    (file, str, file, file, file) -> tuple
    Returns a tuple with the lists of UTR length of miRNA target genes from the TargetScan_mirna file
    that in CNV loss, in CNV gain or not in CNV
    '''

    # get the UTR sequence for each transcripts
    UTRs = UTR_sequence(UTR_file, species)

##    # get all genes in CNV gain and CNV loss
##    CNV_gain, CNV_loss = human_CNV_type(CNV_gene_file, all_gene_file)
    
    # make sets of miRNA transcript targets in CNV gain, CNV loss and non-CNV from CNV_miRNA_file
    transcript_CNV_gain, transcript_CNV_loss, transcript_non_CNV = human_targets_CNV_type(CNV_miRNA_file, CNV_gene_file, all_gene_file)

    # make lists of UTR length for targets in CNV gain, CNV loss and non-CNV
    UTR_CNV_gain = []
    UTR_CNV_loss = []
    UTR_non_CNV = []

    for transcript in UTRs:
        if transcript in transcript_CNV_gain:
            UTR_CNV_gain.append(len(UTRs[transcript]))
        elif transcript in transcript_CNV_loss:
            UTR_CNV_loss.append(len(UTRs[transcript]))
        elif transcript in transcript_non_CNV:
            UTR_non_CNV.append(len(UTRs[transcript]))

    return UTR_CNV_gain, UTR_CNV_loss, UTR_non_CNV
    

def human_cnv_gain_loss_mirna_regulation(CNV_miRNA_file, CNV_gene_file, all_gene_file):
    '''
    (file, file, file) -> dict
    Extract miRNA regulation info for each transcript in CNV_miRNA_file and add the type of CNV (gain, loss, non) for each transcript
    Returns a dictionnary with transcript as key and a list containing the number of targets, the number of sites, the ratio of sites per mirna and the type of CNV
    '''

    # make sets of miRNA transcript targets in CNV gain, CNV loss and non-CNV from CNV_miRNA_file
    transcript_CNV_gain, transcript_CNV_loss, transcript_non_CNV = human_targets_CNV_type(CNV_miRNA_file, CNV_gene_file, all_gene_file)

    # create dictionnary to store mirna regulation information and CNV type status for each transcript
    transcript_CNV_mirna_regulation = {}

    mirnas = open(CNV_miRNA_file, 'r')
    mirnas.readline()
    for line in mirnas:
        if line.rstrip() != '':
            line = line.split()
            transcript = line[1]
            if transcript in transcript_CNV_gain:
                transcript_CNV_mirna_regulation[transcript] = [int(line[2]), int(line[3]), float(line[4]), 'gain']
            elif transcript in transcript_CNV_loss:
                transcript_CNV_mirna_regulation[transcript] = [int(line[2]), int(line[3]), float(line[4]), 'loss']
            elif transcript in transcript_non_CNV:
                transcript_CNV_mirna_regulation[transcript] = [int(line[2]), int(line[3]), float(line[4]), 'non-CNV']

    mirnas.close()
    return transcript_CNV_mirna_regulation


def human_mirna_regulation_CNV_loss_gain_lists(CNV_miRNA_file, CNV_gene_file, all_gene_file):
    '''
    (file, file, file) -> tuple
    Returns a tuple with lists of miRNAS, sites and sites per mirnas for transcripts in file CNV_mirna that are
    in CNV gain, CNV loss or not in CNV
    '''

    # get the miRNA information and CNV type status 
    transcript_CNV_mirna_regulation = human_cnv_gain_loss_mirna_regulation(CNV_miRNA_file, CNV_gene_file, all_gene_file)

    # make lists of variables for gain and loss CNV
    mirnas_gain_CNV = []
    mirnas_loss_CNV = []
    mirnas_non_CNV = []
    sites_gain_CNV = []
    sites_loss_CNV = []
    sites_non_CNV = []
    ratio_gain_CNV = []
    ratio_loss_CNV = []
    ratio_non_CNV = []

    # populate lists
    for transcript in transcript_CNV_mirna_regulation:
        if transcript_CNV_mirna_regulation[transcript][-1] == 'gain':
            mirnas_gain_CNV.append(transcript_CNV_mirna_regulation[transcript][0])
            sites_gain_CNV.append(transcript_CNV_mirna_regulation[transcript][1])
            ratio_gain_CNV.append(transcript_CNV_mirna_regulation[transcript][2])
        elif transcript_CNV_mirna_regulation[transcript][-1] == 'loss':
            mirnas_loss_CNV.append(transcript_CNV_mirna_regulation[transcript][0])
            sites_loss_CNV.append(transcript_CNV_mirna_regulation[transcript][1])
            ratio_loss_CNV.append(transcript_CNV_mirna_regulation[transcript][2])
        elif transcript_CNV_mirna_regulation[transcript][-1] == 'non-CNV':
            mirnas_non_CNV.append(transcript_CNV_mirna_regulation[transcript][0])
            sites_non_CNV.append(transcript_CNV_mirna_regulation[transcript][1])
            ratio_non_CNV.append(transcript_CNV_mirna_regulation[transcript][2])

    return (mirnas_gain_CNV, mirnas_loss_CNV, mirnas_non_CNV,
            sites_gain_CNV, sites_loss_CNV, sites_non_CNV,
            ratio_gain_CNV, ratio_loss_CNV, ratio_non_CNV)


def test_UTR_length_CNV_gain_loss(UTR_file, species, CNV_miRNA_file, CNV_gene_file):
    '''
    (list, list) -> tuple
    Performs a Wilcoxon rank sum test of the UTR length difference between miRNA targets in CNV gain, CNV loss and non-CNV
    Returns a tuple with the z-value and the p-value
    '''
    # get the lists of UTR lengths for genes in CNV gains, CNV losses and not in CNVs
    if species == 'worm':
        gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR = worm_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_gene_file)
    elif species == 'human':
        all_gene_file = 'Homo_sapiens.gene_info'
        gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR = human_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_gene_file, all_gene_file)
    elif species == 'fly':
        CNV_subtype_gene_file = CNV_gene_file
        gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR = fly_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_subtype_gene_file)
 
    # performs the Wilcoxon rank sum test
    from scipy import stats
    diff_UTR_gain_loss = stats.ranksums(gain_CNV_UTR, loss_CNV_UTR)
    diff_UTR_gain_nonCNV = stats.ranksums(gain_CNV_UTR, non_CNV_UTR)
    diff_UTR_loss_nonCNV = stats.ranksums(loss_CNV_UTR, non_CNV_UTR)

    return diff_UTR_gain_loss, diff_UTR_gain_nonCNV, diff_UTR_loss_nonCNV


def average_UTR_length_CNV_gain_loss(UTR_file, species, CNV_miRNA_file, CNV_gene_file):
    '''
    (list, list) -> tuple
    Returns a tuple of 2-item tuples with the mean and standard error of the number of mirnas,
    sites and sites per mirna between CNV gain, CNV loss and non-CNV
    '''
    # get the lists of UTR lengths for genes in CNV gains, CNV losses and not in CNVs
    if species == 'worm':
        gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR = worm_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_gene_file)
    elif species == 'human':
        all_gene_file = 'Homo_sapiens.gene_info'
        gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR = human_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_gene_file, all_gene_file)
    elif species == 'fly':
        CNV_subtype_gene_file = CNV_gene_file
        gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR = fly_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_subtype_gene_file)

    # compute average and SEM
    mean_SEM_UTR_CNV_gain = compute_mean_std_error(gain_CNV_UTR)
    mean_SEM_UTR_CNV_loss = compute_mean_std_error(loss_CNV_UTR)
    mean_SEM_UTR_non_CNV = compute_mean_std_error(non_CNV_UTR)

    return mean_SEM_UTR_CNV_gain, mean_SEM_UTR_CNV_loss, mean_SEM_UTR_non_CNV

    
def test_mirna_regulation_CNV_loss_gain(CNV_miRNA_file, CNV_gene_file, species):
    '''
    (file, file, file) -> tuple
    Returns a tuple of 2-item tuples with the z-score and the p-value comparing the number of mirnas,
    sites and sites per mirna between CNV gain, CNV loss and non-CNV
    '''

    # get the list of mirnas, sites and ratio for CNV gain , CNV loss and non-CNV genes
    if species == 'worm':
        UTR_file = 'Celegans_UTR_Sequences.txt'
        regulation_lists = worm_mirna_regulation_CNV_loss_gain_lists(UTR_file, CNV_miRNA_file, CNV_gene_file)
    elif species == 'human':
        all_gene_file = 'Homo_sapiens.gene_info'
        regulation_lists = human_mirna_regulation_CNV_loss_gain_lists(CNV_miRNA_file, CNV_gene_file, all_gene_file)
    elif species == 'fly':
        CNV_subtype_gene_file = CNV_gene_file 
        regulation_lists = fly_mirna_regulation_CNV_loss_gain_lists(CNV_miRNA_file, CNV_subtype_gene_file)

    mirnas_gain_CNV = regulation_lists[0]
    mirnas_loss_CNV = regulation_lists[1]
    mirnas_non_CNV = regulation_lists[2]
    sites_gain_CNV = regulation_lists[3]
    sites_loss_CNV = regulation_lists[4]
    sites_non_CNV = regulation_lists[5]
    ratio_gain_CNV = regulation_lists[6]
    ratio_loss_CNV = regulation_lists[7]
    ratio_non_CNV = regulation_lists[8]

    # performs Wilcoxon rank sum tests
    from scipy import stats
    diff_mirnas_gain_loss = stats.ranksums(mirnas_gain_CNV, mirnas_loss_CNV)
    diff_sites_gain_loss = stats.ranksums(sites_gain_CNV, sites_loss_CNV)
    diff_ratio_gain_loss = stats.ranksums(ratio_gain_CNV, ratio_loss_CNV)

    diff_mirnas_gain_nonCNV = stats.ranksums(mirnas_gain_CNV, mirnas_non_CNV)
    diff_sites_gain_nonCNV = stats.ranksums(sites_gain_CNV, sites_non_CNV)
    diff_ratio_gain_nonCNV = stats.ranksums(ratio_gain_CNV, ratio_non_CNV)

    diff_mirnas_loss_nonCNV = stats.ranksums(mirnas_loss_CNV, mirnas_non_CNV)
    diff_sites_loss_nonCNV = stats.ranksums(sites_loss_CNV, sites_non_CNV)
    diff_ratio_loss_nonCNV = stats.ranksums(ratio_loss_CNV, ratio_non_CNV)

    return (diff_mirnas_gain_loss, diff_sites_gain_loss, diff_ratio_gain_loss,
            diff_mirnas_gain_nonCNV, diff_sites_gain_nonCNV, diff_ratio_gain_nonCNV,
            diff_mirnas_loss_nonCNV, diff_sites_loss_nonCNV, diff_ratio_loss_nonCNV)


def average_mirna_regulation_CNV_loss_gain(CNV_miRNA_file, CNV_gene_file, species):
    '''
    (file, file, file) -> tuple
    Returns a tuple of 2-item tuples with the mean and standard error of the number of mirnas,
    sites and sites per mirna between CNV gain, CNV loss and non-CNV
    '''

    # get the list of mirnas, sites and ratio for CNV gain , CNV loss and non-CNV genes
    if species == 'worm':
        UTR_file = 'Celegans_UTR_Sequences.txt'
        regulation_lists = worm_mirna_regulation_CNV_loss_gain_lists(UTR_file, CNV_miRNA_file, CNV_gene_file)
    elif species == 'human':
        all_gene_file = 'Homo_sapiens.gene_info'
        regulation_lists = human_mirna_regulation_CNV_loss_gain_lists(CNV_miRNA_file, CNV_gene_file, all_gene_file)
    elif species == 'fly':
        CNV_subtype_gene_file = CNV_gene_file
        regulation_lists = fly_mirna_regulation_CNV_loss_gain_lists(CNV_miRNA_file, CNV_subtype_gene_file)

    mirnas_gain_CNV = regulation_lists[0]
    mirnas_loss_CNV = regulation_lists[1]
    mirnas_non_CNV = regulation_lists[2]
    sites_gain_CNV = regulation_lists[3]
    sites_loss_CNV = regulation_lists[4]
    sites_non_CNV = regulation_lists[5]
    ratio_gain_CNV = regulation_lists[6]
    ratio_loss_CNV = regulation_lists[7]
    ratio_non_CNV = regulation_lists[8]

    # compute average and SEM
    mean_SEM_mirnas_CNV_gain = compute_mean_std_error(mirnas_gain_CNV)
    mean_SEM_mirnas_CNV_loss = compute_mean_std_error(mirnas_loss_CNV)
    mean_SEM_mirnas_non_CNV = compute_mean_std_error(mirnas_non_CNV)
    mean_SEM_sites_CNV_gain = compute_mean_std_error(sites_gain_CNV)
    mean_SEM_sites_CNV_loss = compute_mean_std_error(sites_loss_CNV)
    mean_SEM_sites_non_CNV = compute_mean_std_error(sites_non_CNV)
    mean_SEM_ratio_CNV_gain = compute_mean_std_error(ratio_gain_CNV)
    mean_SEM_ratio_CNV_loss = compute_mean_std_error(ratio_loss_CNV)
    mean_SEM_ratio_non_CNV = compute_mean_std_error(ratio_non_CNV)
   
    return (mean_SEM_mirnas_CNV_gain, mean_SEM_mirnas_CNV_loss, mean_SEM_mirnas_non_CNV,
            mean_SEM_sites_CNV_gain, mean_SEM_sites_CNV_loss, mean_SEM_sites_non_CNV,
            mean_SEM_ratio_CNV_gain, mean_SEM_ratio_CNV_loss, mean_SEM_ratio_non_CNV)


def fly_find_CNV_type(CNV_gene_file, fly_gene_coordinates, CNV_type_file):
    '''
    (file, file, file) -> file
    Save to file the gene and the corresponding type of CNV
    Genes that fall in different types of CNVs (duplication and deletion) have been removed
    '''

    # make a dictionnary FlyBase gene ID: gene name from the CNV_gene_file
    FlyBase_genes = {}
    cnvs = open(CNV_gene_file, 'r')
    cnvs.readline()
    for line in cnvs:
        if line.rstrip() != '':
            line = line.rstrip().split()
            FlyBase_genes[line[0]] = line[1]
    cnvs.close()
    
    # make a dictionnary of gene position for all genes in the genome
    # using the FlyBase gene ID as key, and a list with start, end and chromosome
    genome = {}
    flybase = open(fly_gene_coordinates, 'r')
    for line in flybase:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'FlyBase' and line[2] == 'gene':
                gene = line[8][line[8].index('FBgn'): line[8].index(';')]
                genome[gene] = [line[0], int(line[3]), int(line[4])]
    flybase.close()
    
    # cross reference the dictionnary of gene names and gene positions and 
    # make a dictionnary of gene position for all the CNV genes
    # using the gene name as key, and a list woth start, end and chromosome

    CNV_genes_positions = {}
    for gene in FlyBase_genes:
        gene_name = FlyBase_genes[gene]
        CNV_genes_positions[gene_name] = genome[gene]

    # make a dictionnary with CNV positions and CNV type from the CNV_type_file {i: [start, end, chromosome, type]}
    i = 0
    CNV_type = {}
    cnvs = open(CNV_type_file, 'r')
    cnvs.readline()
    for line in cnvs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            CNV_type[i] = [int(line[0]), int(line[1]), line[2], line[3]]
            i += 1
    cnvs.close()

    # cross reference the dictionnary of CNV gene positions and CNV type
    # make a dictionnary with gene name : CNV type {'CG10345' : 'del'}
    CNV_gene_subtype = {}
    for gene in CNV_genes_positions:
        for i in CNV_type:
            if CNV_genes_positions[gene][0] == CNV_type[i][2]:
                coord1 = set(j for j in range(CNV_genes_positions[gene][1], CNV_genes_positions[gene][2]+1))
                coord2 = set(j for j in range(CNV_type[i][0], CNV_type[i][1]+1))
                if len(coord1.intersection(coord2)) != 0:
                    if gene in CNV_gene_subtype:
                        CNV_gene_subtype[gene].add(CNV_type[i][-1])
                    else:
                        CNV_gene_subtype[gene] = set()
                        CNV_gene_subtype[gene].add(CNV_type[i][-1])
    print('CNV subtype search done')
    print(len(CNV_gene_subtype))


    # remove genes located in different types of CNV
    to_remove = []
    for gene in CNV_gene_subtype:
        if len(CNV_gene_subtype[gene]) > 1:
            to_remove.append(gene)
    for gene in to_remove:
        del CNV_gene_subtype[gene]
    print(len(CNV_gene_subtype))

    # change the type of the value in dict CNV_gene_subtype from set to string
    for gene in CNV_gene_subtype:
        CNV_gene_subtype[gene] = CNV_gene_subtype[gene].pop()

    # open file for saving content of the CNV_gene_subtype dictionnary
    newfile = ('Drosophila_gene_name_CNV_subtype.txt', 'w')
    newfile.write('Gene_name\tCNV_type\n')
    for gene in CNV_gene_subtype:
        newfile.write(gene + '\t' + CNV_gene_subtype[gene] + '\n')
    newfile.close()
    
       

def fly_transcript_to_gene(CNV_miRNA_file):
    '''
    (file) -> dict
    Returns a dictionnary of miRNA target transcript name : gene name pairs extracted from the CNV_miRNA_file
    '''

    # make a dictionnary of transcript : gene name pairs from the CNV_miRNA_file
    transcript_to_gene = {}
    cnvs = open(CNV_miRNA_file, 'r')
    cnvs.readline()
    for line in cnvs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            transcript_to_gene[line[1]] = line[0]
    cnvs.close()
    return transcript_to_gene




def fly_CNV_nonCNV_transcripts(CNV_miRNA_file):
    '''
    (file) -> (set, set)
    Returns a tuple containing the set of miRNA target genes that are in CNVs and the set of genes that are not in CNV
    '''

    # create sets of CNVs and non-CNVs transcripts
    CNVs = set()
    nonCNVs = set()
    
    # open CNV miRNA file for reading
    cnvs = open(CNV_miRNA_file, 'r')
    cnvs.readline()
    for line in cnvs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[-1] == 'CNV':
                CNVs.add(line[1])
            elif line[-1] == 'non-CNV':
                nonCNVs.add(line[1])
    cnvs.close()
    return CNVs, nonCNVs


def fly_grab_CNV_type(CNV_subtype_gene_file):
    '''
    (file) -> dict
    Returns a dictionnary of gene name : CNV type pairs from the the CNV_subtype_gene_file
    Genes that fall in different CNV types (duplication and deletion) have been removed
    '''
    # create a dictionnary to stote the gene name : CNV type pairs
    CNV_types = {}

    # open file for reading
    cnvs = open(CNV_subtype_gene_file, 'r')
    cnvs.readline()
    for line in cnvs:
        line  = line.rstrip()
        if line != '':
            line  =  line.split()
            CNV_types[line[0]] = line[1]
    cnvs.close()
    return CNV_types



def fly_cnv_gain_loss_UTR_length(UTR_file, species, CNV_miRNA_file, CNV_subtype_gene_file):
    '''
    (file, str, file, file) -> tuple
    Returns a tuple with the lists of UTR length of miRNA target genes from the TargetScan_mirna file
    that in CNV loss, in CNV gain or not in CNV
    '''

    # get the UTR sequence of each transcript
    UTR_seq = UTR_sequence(UTR_file, species)

    # make a dictionnary of transcript : gene name pairs from the CNV_miRNA_file
    transcript_to_gene = fly_transcript_to_gene(CNV_miRNA_file)

    # make sets of CNV and non-CNV target transcripts
    CNVs, nonCNVs = fly_CNV_nonCNV_transcripts(CNV_miRNA_file)

    # get the CNV type information for each gene
    CNV_types = fly_grab_CNV_type(CNV_subtype_gene_file)
    

    # store the length of UTR for gain CNVs, loss CNVs and non-CNVs
    gain_CNV_UTR = []
    loss_CNV_UTR = []
    non_CNV_UTR = []
    for transcript in transcript_to_gene:
        UTR_length = len(UTR_seq[transcript])
        gene = transcript_to_gene[transcript]
        # genes in different CNV subtypes have been removed
        if transcript in CNVs and gene in CNV_types:
            if CNV_types[gene] == 'dup':
                gain_CNV_UTR.append(UTR_length)
            elif CNV_types[gene] == 'del':
                loss_CNV_UTR.append(UTR_length)
        elif transcript in nonCNVs:
            non_CNV_UTR.append(UTR_length)


    return gain_CNV_UTR, loss_CNV_UTR, non_CNV_UTR


def fly_mirnas_cnv_gain_loss(CNV_miRNA_file, CNV_subtype_gene_file):
    '''
    (file, str, file, file) -> dict
    Returns a dictionnary with transcript as key and a list with the number of mirnas,
    number of sites, ratio and CNV type status from the CNV_miRNA_file
    '''

    # get the CNV type information for each gene
    CNV_types = fly_grab_CNV_type(CNV_subtype_gene_file)

    # make a dictionnary of mirna regulation with gene as key
    cnv_regulation = {}
    cnv_mirnas = open(CNV_miRNA_file, 'r')
    cnv_mirnas.readline()
    for line in cnv_mirnas:
        if line.rstrip() != '':
            line = line.split()
            gene = line[0]
            if gene in CNV_types and line[-1] == 'CNV':
                if CNV_types[gene] == 'dup':
                    cnv_regulation[gene] = [int(line[2]), int(line[3]), float(line[4]), 'dup']
                elif CNV_types[gene] == 'del':
                    cnv_regulation[gene] = [int(line[2]), int(line[3]), float(line[4]), 'del']
            elif line[-1] == 'non-CNV':
                cnv_regulation[gene] = [int(line[2]), int(line[3]), float(line[4]), 'non-CNV']
                
    cnv_mirnas.close()
    return cnv_regulation


def fly_mirna_regulation_CNV_loss_gain_lists(CNV_miRNA_file, CNV_subtype_gene_file):
    '''
    (file, file, file) -> tuple
    Returns a tuple with lists of miRNAS, sites and sites per mirnas for transcripts in file CNV_mirna that are
    in CNV gain, CNV loss or not in CNV
    '''

    # get the miRNA information and CNV type status 
    cnv_regulation = fly_mirnas_cnv_gain_loss(CNV_miRNA_file, CNV_subtype_gene_file)

    # make lists of variables for gain and loss CNV
    mirnas_gain_CNV = []
    mirnas_loss_CNV = []
    mirnas_non_CNV = []
    sites_gain_CNV = []
    sites_loss_CNV = []
    sites_non_CNV = []
    ratio_gain_CNV = []
    ratio_loss_CNV = []
    ratio_non_CNV = []

    # populate lists
    for gene in cnv_regulation:
        if cnv_regulation[gene][-1] == 'dup':
            mirnas_gain_CNV.append(cnv_regulation[gene][0])
            sites_gain_CNV.append(cnv_regulation[gene][1])
            ratio_gain_CNV.append(cnv_regulation[gene][2])
        elif cnv_regulation[gene][-1] == 'del':
            mirnas_loss_CNV.append(cnv_regulation[gene][0])
            sites_loss_CNV.append(cnv_regulation[gene][1])
            ratio_loss_CNV.append(cnv_regulation[gene][2])
        elif cnv_regulation[gene][-1] == 'non-CNV':
            mirnas_non_CNV.append(cnv_regulation[gene][0])
            sites_non_CNV.append(cnv_regulation[gene][1])
            ratio_non_CNV.append(cnv_regulation[gene][2])

    return (mirnas_gain_CNV, mirnas_loss_CNV, mirnas_non_CNV,
            sites_gain_CNV, sites_loss_CNV, sites_non_CNV,
            ratio_gain_CNV, ratio_loss_CNV, ratio_non_CNV)
