from CNV_miRNAs_TargetScan import *


def compute_quartiles(L):
    '''
    (list) -> tuple
    Return the three quartile points of values contained in list L
    
    >>> compute_quartiles([1, 3, 3, 4, 5, 6, 6, 7, 8, 8,])
    (3, 5.5, 7)
    >>> compute_quartiles([3, 4, 4, 5, 6, 8, 8])
    (4, 5, 8)
    '''
  
    L.sort()
    if len(L) % 2 == 1:
        pos_median = int(len(L) / 2)
        median = L[pos_median]
        first_half = L[:pos_median]
        second_half = L[pos_median + 1:]

    elif len(L) % 2 == 0:
        pos1 = int(len(L) / 2)
        pos2 = pos1 - 1
        median = (L[pos1] + L[pos2]) / 2
        first_half = L[:pos1]
        second_half = L[pos1:]

    pos_Q1 = int(len(first_half) / 2)
    Q1 = first_half[pos_Q1]
    pos_Q3 = int(len(second_half) / 2)
    Q3 = second_half[pos_Q3] 

    return (Q1, median, Q3)


def human_CNV_miRNA_table_sorted_UTR_length(UTR_file, species, summary_counts, CNV_file, all_gene_file,
                                            miRNA_family_info, conservation_sites):
    '''
    (file, str, file, file, file, file, str) -> file, file, file, file
    For a single transcript / gene sorted by 3' UTR length quartile , write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not. Write 4 files for each length quartiles
    '''

    # get the UTR sequence for each transcript
    transcripts_UTR = UTR_sequence(UTR_file, species)
    print('get UTR done')

    # get the gene ID for each target transcript
    transcripts = transcript_gene_pairs(summary_counts)
    print('get transcript name done')

    # get CNV status
    CNV_genes = human_CNV_genes(CNV_file, all_gene_file)
    print('get CNV status done')

    # get valid genes
    valid_genes = sort_valid_human_genes(all_gene_file)
    print('get valid genes done')

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)
    print('get mirna family info done')

    # get mirna regularion for each transcript
    regulators = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)
    print('get mirna regulation info for each transcript done')

    # create a dictionnary to store info about miRNA regultors for each gene, transcript_ID, length of UTR and CNV status
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, length_UTR, CNV_status]}
    gene_regulated = {}
    for transcript_ID in regulators:
        gene = transcripts[transcript_ID]
        UTR_length = len(transcripts_UTR[transcript_ID])
        if gene not in gene_regulated:
            gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
            gene_regulated[gene].append(transcript_ID) # add the transcript ID
            gene_regulated[gene].append(UTR_length) # get the UTR length of the given transcript

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

    print('get mirna info for each done done')

    # make a list of UTR lenth recorded in gene_regulated
    all_UTR = []
    for gene in gene_regulated:
        all_UTR.append(gene_regulated[gene][-2])

    # compute length quartiles
    length_quartiles = compute_quartiles(all_UTR)

    # create dictionnaries of with mirna regulation infor for each gene sorted according to the length of their transcript
    Q1, Q2, Q3, Q4 = {}, {}, {}, {}

    # sort genes according to UTR length
    for gene in gene_regulated:
        UTR_length = gene_regulated[gene][-2]
        if UTR_length < length_quartiles[0]:
            Q1[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[0] and UTR_length < length_quartiles[1]:
            Q2[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[1] and UTR_length < length_quartiles[2]:
            Q3[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[2]:
            Q4[gene] = gene_regulated[gene]
    print('UTR length sorted in length quartiles done')

    # store the dictionnaries into a list to loop over
    Q_length = [Q1, Q2, Q3, Q4]
    filenames = ['Human_miRNAs_TargetScan_allsites_Q1.txt', 'Human_miRNAs_TargetScan_allsites_Q2.txt',
                 'Human_miRNAs_TargetScan_allsites_Q3.txt', 'Human_miRNAs_TargetScan_allsites_Q4.txt']

    # go through each quartile dictionnary, and write in corresponding file
    for i in range(len(Q_length)):
        # write to file
        newfile = open(filenames[i], 'w')
        newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
        for gene in Q_length[i]:
            newfile.write(gene + '\t') # write gene name
            newfile.write(Q_length[i][gene][-3] + '\t') # write transcript name
            N_mirnas = 0
            N_sites = 0
            for pair in Q_length[i][gene][:-3]: # compute number of mirnas and number of sites 
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/ N_mirnas) + '\t' + Q_length[i][gene][-1] + '\n')
        newfile.close()


def worm_CNV_miRNA_table_sorted_UTR_length(UTR_file, species, summary_counts, CNV_file, all_gene_file,
                                           miRNA_family_info, conservation_sites):
    '''
    (file, str, file, file, file, file, str) -> file, file, file, file
    For a single transcript / gene sorted by 3' UTR length quartile , write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not. Write 4 files for each length quartiles
    '''

    # get the UTR sequence for each transcript
    transcripts_UTR = UTR_sequence(UTR_file, species)
    print('get UTR done')

    # get the gene ID for each target transcript
    transcripts = transcript_gene_pairs(summary_counts)
    print('get transcript name done')

    # get CNV status
    CNV_genes = worm_CNV_genes(CNV_file, all_gene_file)
    print('get CNV status done')

    # get valid genes
    valid_genes = sort_valid_worm_genes(all_gene_file)
    print('get valid genes done')

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)
    print('get mirna family info done')

    # get mirna regularion for each transcript
    regulators = worm_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)
    print('get mirna regulation info for each transcript done')

    # create a dictionnary to store info about miRNA regultors for each gene, transcript_ID, length of UTR and CNV status
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, length_UTR, CNV_status]}
    gene_regulated = {}
    for transcript_ID in regulators:
        gene = transcripts[transcript_ID]
        UTR_length = len(transcripts_UTR[transcript_ID])
        if gene not in gene_regulated:
            gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
            gene_regulated[gene].append(transcript_ID) # add the transcript ID
            gene_regulated[gene].append(UTR_length) # get the UTR length of the given transcript

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

    print('get mirna info for each done done')

    # make a list of UTR lenth recorded in gene_regulated
    all_UTR = []
    for gene in gene_regulated:
        all_UTR.append(gene_regulated[gene][-2])

    # compute length quartiles
    length_quartiles = compute_quartiles(all_UTR)

    # create dictionnaries of with mirna regulation infor for each gene sorted according to the length of their transcript
    Q1, Q2, Q3, Q4 = {}, {}, {}, {}

    # sort genes according to UTR length
    for gene in gene_regulated:
        UTR_length = gene_regulated[gene][-2]
        if UTR_length < length_quartiles[0]:
            Q1[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[0] and UTR_length < length_quartiles[1]:
            Q2[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[1] and UTR_length < length_quartiles[2]:
            Q3[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[2]:
            Q4[gene] = gene_regulated[gene]
    print('UTR length sorted in length quartiles done')

    # store the dictionnaries into a list to loop over
    Q_length = [Q1, Q2, Q3, Q4]
    filenames = ['Worm_miRNAs_TargetScan_allsites_Q1.txt', 'Worm_miRNAs_TargetScan_allsites_Q2.txt',
                 'Worm_miRNAs_TargetScan_allsites_Q3.txt', 'Worm_miRNAs_TargetScan_allsites_Q4.txt']

    # go through each quartile dictionnary, and write in corresponding file
    for i in range(len(Q_length)):
        # write to file
        newfile = open(filenames[i], 'w')
        newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
        for gene in Q_length[i]:
            newfile.write(gene + '\t') # write gene name
            newfile.write(Q_length[i][gene][-3] + '\t') # write transcript name
            N_mirnas = 0
            N_sites = 0
            for pair in Q_length[i][gene][:-3]: # compute number of mirnas and number of sites 
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/ N_mirnas) + '\t' + Q_length[i][gene][-1] + '\n')
        newfile.close()


def fly_CNV_miRNA_table_sorted_UTR_length(UTR_file, species, CNV_genes_file, fly_gene_coordinates,
                                          miRNA_family_info, summary_counts, conservation_sites):
    '''
    (file, str, file, file, file, file, str) -> file, file, file, file
    For a single transcript / gene sorted by 3' UTR length quartile , write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not. Write 4 files for each length quartiles
    '''

    # get the UTR sequence for each transcript
    transcripts_UTR = UTR_sequence(UTR_file, species)
    print('get UTR done')

    # get the gene ID for each target transcript
    transcripts = fly_transcript_gene_pairs(UTR_file)
    print('get transcript name done')

    # get the CNV genes
    CNV_genes = set()
    CNV = open(CNV_genes_file, 'r')
    CNV.readline()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            line = line.split()
            CNV_genes.add(line[0])
            CNV_genes.add(line[1])
    CNV.close()
    print('get CNV status done')

    # get valid genes
    valid_genes = sort_valid_fly_genes(fly_gene_coordinates)
    print('get valid genes done')

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)
    print('get mirna family info done')

    # get mirna regularion for each transcript
    regulators = fly_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)
    print('get mirna regulation info for each transcript done')

    # create a dictionnary to store info about miRNA regultors for each gene, transcript_ID, length of UTR and CNV status
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, length_UTR, CNV_status]}
    gene_regulated = {}
    for transcript_ID in regulators:
        gene = transcripts[transcript_ID]
        UTR_length = len(transcripts_UTR[transcript_ID])
        if gene not in gene_regulated:
            gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
            gene_regulated[gene].append(transcript_ID) # add the transcript ID
            gene_regulated[gene].append(UTR_length) # get the UTR length of the given transcript

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

    print('get mirna info for each done done')

    # make a list of UTR lenth recorded in gene_regulated
    all_UTR = []
    for gene in gene_regulated:
        all_UTR.append(gene_regulated[gene][-2])

    # compute length quartiles
    length_quartiles = compute_quartiles(all_UTR)

    # create dictionnaries of with mirna regulation infor for each gene sorted according to the length of their transcript
    Q1, Q2, Q3, Q4 = {}, {}, {}, {}

    # sort genes according to UTR length
    for gene in gene_regulated:
        UTR_length = gene_regulated[gene][-2]
        if UTR_length < length_quartiles[0]:
            Q1[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[0] and UTR_length < length_quartiles[1]:
            Q2[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[1] and UTR_length < length_quartiles[2]:
            Q3[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[2]:
            Q4[gene] = gene_regulated[gene]
    print('UTR length sorted in length quartiles done')

    # store the dictionnaries into a list to loop over
    Q_length = [Q1, Q2, Q3, Q4]
    filenames = ['Fly_miRNAs_TargetScan_allsites_Q1.txt', 'Fly_miRNAs_TargetScan_allsites_Q2.txt',
                 'Fly_miRNAs_TargetScan_allsites_Q3.txt', 'Fly_miRNAs_TargetScan_allsites_Q4.txt']

    # go through each quartile dictionnary, and write in corresponding file
    for i in range(len(Q_length)):
        # write to file
        newfile = open(filenames[i], 'w')
        newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
        for gene in Q_length[i]:
            newfile.write(gene + '\t') # write gene name
            newfile.write(Q_length[i][gene][-3] + '\t') # write transcript name
            N_mirnas = 0
            N_sites = 0
            for pair in Q_length[i][gene][:-3]: # compute number of mirnas and number of sites 
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/ N_mirnas) + '\t' + Q_length[i][gene][-1] + '\n')
        newfile.close()




def zebrafish_CNV_miRNA_table_sorted_UTR_length(UTR_file, species, CNV_genes_file, miRNA_family_info, summary_counts,
                                                ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only = 'yes'):
    '''
    (file, str, file, file, file, file, str) -> file, file, file, file
    For a single transcript / gene sorted by 3' UTR length quartile , write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not. Write 4 files for each length quartiles
    '''

    # get the UTR sequence for each transcript
    transcripts_UTR = UTR_sequence(UTR_file, species)
    print('get UTR done')

    # get the CNV genes
    CNV = open(CNV_genes_file, 'r')
    CNV_genes = set()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            CNV_genes.add(line)
    CNV.close()
    print('get CNV status done')

    # get valid genes
    valid_genes = sort_valid_zebrafish_genes(ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only)
    print('get valid genes done')

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)
    print('get mirna family info done')

    # get mirna regularion for each transcript
    regulators = zebrafish_miRNA_regulation_TargetScan_sites(summary_counts, species, miRNA_family_info)
    print('get mirna regulation info for each transcript done')

    # create a dictionnary to store info about miRNA regultors for each gene, transcript_ID, length of UTR and CNV status
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, length_UTR, CNV_status]}
    gene_regulated = {}
    for transcript_ID in regulators:
        gene = transcript_ID[:transcript_ID.index('.')] # the transcript ID is not a valid ensembl ID but the gene_ID + .number
        UTR_length = len(transcripts_UTR[transcript_ID])
        if gene not in gene_regulated:
            gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
            gene_regulated[gene].append(transcript_ID) # add the transcript ID
            gene_regulated[gene].append(UTR_length) # get the UTR length of the given transcript

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

    print('get mirna info for each done done')

    # make a list of UTR lenth recorded in gene_regulated
    all_UTR = []
    for gene in gene_regulated:
        all_UTR.append(gene_regulated[gene][-2])

    # compute length quartiles
    length_quartiles = compute_quartiles(all_UTR)

    # create dictionnaries of with mirna regulation infor for each gene sorted according to the length of their transcript
    Q1, Q2, Q3, Q4 = {}, {}, {}, {}

    # sort genes according to UTR length
    for gene in gene_regulated:
        UTR_length = gene_regulated[gene][-2]
        if UTR_length < length_quartiles[0]:
            Q1[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[0] and UTR_length < length_quartiles[1]:
            Q2[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[1] and UTR_length < length_quartiles[2]:
            Q3[gene] = gene_regulated[gene]
        elif UTR_length >= length_quartiles[2]:
            Q4[gene] = gene_regulated[gene]
    print('UTR length sorted in length quartiles done')

    # store the dictionnaries into a list to loop over
    Q_length = [Q1, Q2, Q3, Q4]
    filenames = ['Zebrafish_miRNAs_TargetScan_allsites_Q1.txt', 'Zebrafish_miRNAs_TargetScan_allsites_Q2.txt',
                 'Zebrafish_miRNAs_TargetScan_allsites_Q3.txt', 'Zebrafish_miRNAs_TargetScan_allsites_Q4.txt']

    # go through each quartile dictionnary, and write in corresponding file
    for i in range(len(Q_length)):
        # write to file
        newfile = open(filenames[i], 'w')
        newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
        for gene in Q_length[i]:
            newfile.write(gene + '\t') # write gene name
            newfile.write(Q_length[i][gene][-3] + '\t') # write transcript name
            N_mirnas = 0
            N_sites = 0
            for pair in Q_length[i][gene][:-3]: # compute number of mirnas and number of sites 
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/ N_mirnas) + '\t' + Q_length[i][gene][-1] + '\n')
        newfile.close()

