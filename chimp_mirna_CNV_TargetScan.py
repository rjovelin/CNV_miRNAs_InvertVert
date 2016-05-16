from CNV_miRNAs_TargetScan import *


def chimp_UTR_sequence(UTR_file):
    '''
    (file, str) -> dict
    Return a dictionnary with the UTR sequence of each transcript
    '''

    # open file for reading
    utr = open(UTR_file, 'r')
    header = utr.readline()

    # create a dictionnary to store the urt sequence for each transcript
    utr_sequences = {}

    # read the file
    for line in utr:
        line = line.rstrip()
        if line != '':
            line = line.split()
            transcript = line[0]
            sequence = line[-1]
            while '-' in sequence:
                sequence = sequence.replace('-', '')
            if '9598' in line:
                utr_sequences[transcript] = sequence

    utr.close()
    return utr_sequences


def chimp_miRNA_family(miRNA_family_info):
    '''
    (file, str) -> dict
    Return a dictionnary with seed as key and a list of same-family miRNAs as value
    '''

    family = {}
    miRNA = open(miRNA_family_info, 'r')
    header = miRNA.readline()

    for line in miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[2] == '9598':
                if line[1] in family:
                    family[line[1]].append(line[3])
                else:
                    family[line[1]] = [line[3]]
    miRNA.close()
    return family


def chimp_transcript_to_gene(ensemblToGeneName):
    '''
    (file) -> dict
    Return a dictionnary with ensembl transcript ID as key and common gene ID as value
    '''
    
    # make a dictionnary with the ensembl transcript : gene pairs
    transcripts_genes = {}
    transcripts = open(ensemblToGeneName, 'r')
    transcripts.readline()
    for line in transcripts:
        line = line.rstrip()
        if line !='':
            line = line.split()
            transcripts_genes[line[0]] = line[1]

    transcripts.close()
    return transcripts_genes


def chimp_find_CNV_genes(CNV_positions_file, ensembl_transcript_coordinates, ensemblToGeneName):
    '''
    (file, file, file) -> file
    Save to file the set of chimp CNV genes
    '''

    # make a dictionnary to store the coordinates of each CNV {CNV1:[chromo, start, end]}
    cnv_coord = open(CNV_positions_file, 'r')
    cnv_coord.readline()
    cnv_coord.readline()
    CNV_positions = {}
    i = 0
    for line in cnv_coord:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[0].startswith('Chimp'):
                CNV_positions[i] = ['chr' + line[1], int(line[2]), int(line[3])]
                i += 1

    # make a dictionnary with the ensembl transcript : gene pairs
    transcripts_genes = chimp_transcript_to_gene(ensemblToGeneName)
    
    # make a dictionnary to store the coordinates of the ensembl transcripts {TS1 : [chromo, start, end]}
    transcripts_coord = open(ensembl_transcript_coordinates, 'r')
    transcripts_coord.readline()
    transcripts_positions = {}
    for line in transcripts_coord:
        line = line.rstrip()
        if line !='':
            line = line.split()
            transcripts_positions[line[1]] = [line[2], int(line[4]), int(line[5])]

    print(len(transcripts_positions))
    
    # search for CNV transcripts using the transcript and CNV coordinates and store in a set
    # report transcripts that overlap with a CNV even in the UTR or non-coding regions
    done = 0
    chimp_CNV_genes = set()
    for transcript in transcripts_positions:
        if done % 5 == 0:
            print(done, len(chimp_CNV_genes), sep = '\t')
        done += 1
        if transcript in transcripts_genes:
            if transcripts_genes[transcript] not in chimp_CNV_genes:     # no need to search again if gene already CNV
                for CNV in CNV_positions:
                    if transcripts_positions[transcript][0] == CNV_positions[CNV][0]:
                        ts_pos = set(range(transcripts_positions[transcript][1], transcripts_positions[transcript][2] + 1))
                        cnv_pos = set(range(CNV_positions[CNV][1], CNV_positions[CNV][2] + 1))
                        # if CNV and transcript coordinates overlap keep transcript
                        if len(ts_pos.intersection(cnv_pos)) != 0:
                            chimp_CNV_genes.add(transcripts_genes[transcript]) 
                            break                                            # no need to search again if gene already CNV
        
    cnv_coord.close()
    transcripts_coord.close()
    return chimp_CNV_genes


# declare valid_chimp_chromos as global variable
# valid chromosomes are the chromosomes in the ensemble_transcripts_coordinates
# file that are used to scan for CNV genes
valid_chimp_chromos = set('chr' + str(i) for i in range(1, 23))


def chimp_valid_genes(ensembl_transcript_coordinates, ensemblToGeneName):
    '''
    (file, file) -> set
    Returns a set with all genes in chimp located on the same chromosomes scanned to find CNV genes
    '''

    # make a set of transcripts that are located on valid chromosomes
    transcripts = set()
    transcripts_coord = open(ensembl_transcript_coordinates, 'r')
    transcripts_coord.readline()
    for line in transcripts_coord:
        line = line.rstrip()
        if line !='':
            line = line.split()
            if line[2] in valid_chimp_chromos:
                transcripts.add(line[1])

    # make a dictionnary of transcripts : gene name pairs
    transcripts_genes = chimp_transcript_to_gene(ensemblToGeneName)

    # make a set of valid genes
    valid_genes = set()
    for ts_name in transcripts:
        if ts_name in transcripts_genes:
            valid_genes.add(transcripts_genes[ts_name])

    transcripts_coord.close()
    return valid_genes
    
    
def chimp_miRNA_regulation_TargetScan_sites(summary_counts, conservation_sites):
    '''
    (file, str) -> dict
    Return a dictionnary with the number of miRNA regulators, miRNA binding sites and number of sites per miRNA for each target transcript
    exctracted from the summary_count file. Use either only conserved sites or all sites.
    '''

    # open the summary_count file for reading
    summary = open(summary_counts, 'r')
    header = summary.readline()

    # create a dictionnary to store the miRNA family and number of sites for each transcripts
    # {transcript_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2]]}
    regulated_transcripts = {}

    # read the summary_count file
    for line in summary:
        line = line.rstrip()
        if line != '':
            line = line.split()
            transcript = line[0]
            family = line[2]
            N_conserved_sites = int(line[4])
            N_poorly_conserved_sites = int(line[8])
            if line[3] == '9598':
                if conservation_sites == 'all_sites':
                    if transcript in regulated_transcripts:
                        regulated_transcripts[transcript].append([family, (N_conserved_sites + N_poorly_conserved_sites)])
                    else:
                        regulated_transcripts[transcript] = [[family, (N_conserved_sites + N_poorly_conserved_sites)]]
                elif conservation_sites == 'conserved_sites':
                    if N_conserved_sites > 0:
                        if transcript in regulated_transcripts:
                            regulated_transcripts[transcript].append([family, N_conserved_sites])
                        else:
                            regulated_transcripts[transcript] = [[family, N_conserved_sites]]
                
    summary.close()
    return regulated_transcripts
    

def human_transcripts_targets(human_CNV_miRNA_file):
    '''
    (file) -> (set, set)
    Returns the set of human transcripts and the set of corresponding genes targeted by miRNAs from the file CNV_miRNA
    '''
    # create set of transcripts
    transcripts = set()
    genes = set()
    
    # open the file for reading
    human = open(human_CNV_miRNA_file, 'r')
    human.readline()
    for line in human:
        line = line.rstrip()
        if line !='':
            line = line.split()
            transcripts.add(line[1])
            genes.add(line[0])
    human.close()
    return transcripts, genes
            
def chimp_make_miRNA_regulators_table(summary_counts, conservation_sites, miRNA_family_info, ensembl_transcript_coordinates,
                                      ensemblToGeneName, CNV_genes_file, human_CNV_miRNA_file, outputfile):
    '''
    (file, str, file, file, file, file, file) -> file
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not
    '''

    # get the gene ID for each target transcript
    transcripts = transcript_gene_pairs(summary_counts)

    # get the list of miRNAs for each family
    family = chimp_miRNA_family(miRNA_family_info)

    # get the human transcripts regulated by miRNAs
    human_transcripts, human_genes = human_transcripts_targets(human_CNV_miRNA_file)

    # get the CNV status
    CNV_genes = set()
    cnv_file  = open(CNV_genes_file, 'r')
    for line in cnv_file:
        line = line.rstrip()
        if line != '':
            CNV_genes.add(line)
    cnv_file.close()

    # get the miRNA regulation
    regulators = chimp_miRNA_regulation_TargetScan_sites(summary_counts, conservation_sites)
            
    # create a dictionnary to store info about miRNA regultors for each gene
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, CNV_status]}

    gene_regulated = {}
    for transcript_ID in regulators:
        gene = transcripts[transcript_ID]
        if transcript_ID in human_transcripts: # use the same transcripts as in human for direct comparison
            gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
            gene_regulated[gene].append(transcript_ID) # add the transcript ID
        elif transcript_ID not in human_transcripts:
            if gene not in gene_regulated and gene not in human_genes: # if gene in human_genes then skipped because another transcript should be used
                gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
                gene_regulated[gene].append(transcript_ID) # add the transcript ID

    # add the CNV status
    for gene in gene_regulated:
        if gene in CNV_genes:
            gene_regulated[gene].append('CNV')
        else:
            gene_regulated[gene].append('non-CNV')

    # get the set of valid genes
    valid_genes = chimp_valid_genes(ensembl_transcript_coordinates, ensemblToGeneName)

    # remove genes that are not located on the same chromosomes used to scan CNV genes
    to_remove = []
    for gene in gene_regulated:
        if gene not in valid_genes:
            to_remove.append(gene)
    for gene in to_remove:
        del gene_regulated[gene]
    

    # write to file
    newfile = open(outputfile, 'w')
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



def make_miRNA_regulators_table_conserved_families(focal_species, conservation_sites, human_CNV_miRNA_file, outputfile, summary_counts = 'Human_Summary_Counts.txt',
                                                   miRNA_family_info = 'Human_miR_Family_Info.txt', species = 'human',
                                                   CNV_file = 'GRCh37_hg19_variants_2013-05-31.txt', all_gene_file = 'Homo_sapiens.gene_info',
                                                   ensembl_transcript_coordinates = 'panTro2_ensGene', ensemblToGeneName = 'panTro2_ensemblToGeneName',
                                                   CNV_genes_file = 'Chimp_CNV_genes.txt'):
    '''
    (str, str, file, file, file, file, str, file, file, file, file, file)
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not
    Use only conserved miRNA families between human and chimp
    '''

    # get the gene ID for each target transcript
    transcripts = transcript_gene_pairs(summary_counts)

    # get the list of miRNAs for each family in human and in chimp
    if focal_species == 'chimp':
        family = chimp_miRNA_family(miRNA_family_info)
        other_family = miRNA_family(miRNA_family_info, species)
    elif focal_species == 'human':
        family = miRNA_family(miRNA_family_info, species)
        other_family = chimp_miRNA_family(miRNA_family_info)
    
    # keep only conserved families between chimp and human
    family_to_remove = []
    for seed in family:
        if seed not in other_family:
            family_to_remove.append(seed)
    for seed in family_to_remove:
        del family[seed]

    # get the CNV status and sort the valid genes
    # get the miRNA regulators and binding sites
    if focal_species == 'human':
        CNV_genes = human_CNV_genes(CNV_file, all_gene_file)
        valid_genes = sort_valid_human_genes(all_gene_file)
        regulators = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)
    elif focal_species == 'chimp':
        # get the CNV status
        CNV_genes = set()
        cnv_file  = open(CNV_genes_file, 'r')
        for line in cnv_file:
            line = line.rstrip()
            if line != '':
                CNV_genes.add(line)
        cnv_file.close()
        valid_genes = chimp_valid_genes(ensembl_transcript_coordinates, ensemblToGeneName)
        regulators = chimp_miRNA_regulation_TargetScan_sites(summary_counts, conservation_sites)
        # get the human transcripts regulated by miRNAs
        human_transcripts, human_genes = human_transcripts_targets(human_CNV_miRNA_file)


    # create a dictionnary to store info about miRNA regultors for each gene
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, CNV_status]}
    gene_regulated = {}

    if focal_species == 'human':
        for transcript_ID in regulators:
            gene = transcripts[transcript_ID]
            if gene not in gene_regulated:
                gene_regulated[gene] = regulators[transcript_ID] # get the info for a single transcript
                gene_regulated[gene].append(transcript_ID) # add the transcript ID
    elif focal_species == 'chimp':
        for transcript_ID in regulators:
            gene = transcripts[transcript_ID]
            if transcript_ID in human_transcripts: # use the same transcripts as in human for direct comparison
                gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
                gene_regulated[gene].append(transcript_ID) # add the transcript ID
            elif transcript_ID not in human_transcripts:
                if gene not in gene_regulated and gene not in human_genes: # if gene in human_genes then skipped because another transcript should be used
                    gene_regulated[gene] = list(regulators[transcript_ID]) # get the info for a single transcript
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
  
    # write to file
    newfile = open(outputfile, 'w')
    newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
    for gene in gene_regulated:
        N_mirnas = 0
        N_sites = 0
        for pair in gene_regulated[gene][:-2]:
            if pair[0] in family:
                N_mirnas += len(family[pair[0]])
                N_sites += pair[1]
        if N_mirnas != 0:
            newfile.write(gene + '\t')
            newfile.write(gene_regulated[gene][-2] + '\t')
            newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/ N_mirnas) + '\t' + gene_regulated[gene][-1] + '\n')

    newfile.close()



def get_CNV_genes_Perry_study(CNV_file, all_gene_file):
    '''
    (file, file) -> set
    Returns a set of CNV genes from the Perry study
    '''

    # get valid human genes
    valid_genes = sort_valid_human_genes(all_gene_file)

    # make a set to store the CNV genes
    # Note that the affected genes are sometimes line[-1] when no sample ID is not provided
    # or line[-2] when sample ID is provided: take all and filter out the non valid genes

    CNV = open(CNV_file, 'r')
    CNV_genes = set()
    header = CNV.readline()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            if line[4] == 'CNV' and line[7] == '18775914':
                if ',' in line[-2]:
                    genes = line[-2].split(',')
                    for gene in genes:
                        CNV_genes.add(gene)
                elif ',' not in line[-2]:
                    CNV_genes.add(line[-2])
                if ',' in line[-1]:
                    genes = line[-1].split(',')
                    for gene in genes:
                        CNV_genes.add(gene)
                elif ',' not in line[-1]:
                    CNV_genes.add(line[-1])

    # remove non valid genes
    to_remove = []
    for gene in CNV_genes:
        if gene not in valid_genes:
            to_remove.append(gene)
    for gene in to_remove:
        CNV_genes.discard(gene)

    CNV.close()
    return CNV_genes


def get_human_chimp_orthologs(human_CNV_miRNA_file, chimp_CNV_miRNA_file):
    '''
    (file, file) -> set
    Returns a set of human and chimp orthologs regulated by miRNAs
    '''

    # make a set of orthologous genes
    human_CNV_miRNA = open(human_CNV_miRNA_file, 'r')
    human_genes = set()
    for line in human_CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            human_genes.add(line[0])
    human_CNV_miRNA.close()

    chimp_CNV_miRNA = open(chimp_CNV_miRNA_file, 'r')
    chimp_genes = set()
    for line in chimp_CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            chimp_genes.add(line[0])
    chimp_CNV_miRNA.close()

    orthologs = human_genes.intersection(chimp_genes)

    return orthologs


def get_human_chimp_miRNA_regulation(human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only):
    '''
    (file, file, file, str) -> tuple
    Returns a 12-item tuple containing the lists of the number of miRNA regulators,
    the number of miRNA binding sites, and the number of sites per miRN for CNV and non-CNV orthologous genes in human and chimp
    Option to include all human CNV genes or only human CNV genes from the same study identifying chimp CNVs
    The lists must be ordered to perform a paired test
    '''

    # make a set of orthologous genes
    orthologs = get_human_chimp_orthologs(human_CNV_miRNA_file, chimp_CNV_miRNA_file)
     
    # open files for reading
    human_CNV_miRNA = open(human_CNV_miRNA_file, 'r')
    chimp_CNV_miRNA = open(chimp_CNV_miRNA_file, 'r')

    # make dictionnaries with gene as key and the list of mirna information as value
    human_CNV_miRNA.readline()
    chimp_CNV_miRNA.readline()

    human_mirna_regulation = {}
    for line in human_CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            human_mirna_regulation[line[0]] = [int(line[2]), int(line[3]), float(line[4]), line[-1]]
    human_CNV_miRNA.close()

    chimp_mirna_regulation = {}
    for line in chimp_CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            chimp_mirna_regulation[line[0]] = [int(line[2]), int(line[3]), float(line[4]), line[-1]]
    chimp_CNV_miRNA.close()

    # remove genes that are not orthologs
    human_to_remove = []
    for gene in human_mirna_regulation:
        if gene not in orthologs:
            human_to_remove.append(gene)
    for gene in human_to_remove:
        del human_mirna_regulation[gene]

    chimp_to_remove = []
    for gene in chimp_mirna_regulation:
        if gene not in orthologs:
            chimp_to_remove.append(gene)
    for gene in chimp_to_remove:
        del chimp_mirna_regulation[gene]

    # use only the human CNV genes from the same study identifying chimp CNV genes or use all human CNV genes
    # if perry_only == yes : use only human and chimp CNV genes from the same study
    if perry_only == 'yes':
        perry_human_CNV_genes = get_CNV_genes_Perry_study(CNV_file, all_gene_file)
        # change the CNV status of the human target genes
        for gene in human_mirna_regulation:
            if gene in perry_human_CNV_genes:
                human_mirna_regulation[gene][-1] = 'CNV'
            elif gene not in perry_human_CNV_genes:
                human_mirna_regulation[gene][-1] = 'non-CNV'
        
    # create ordered lists to store the values of the different metrics for CNV and non-CNV genes
    human_CNV_mirnas = [] # store N mirnas for human CNV genes
    human_CNV_sites = []
    human_CNV_ratio = []
    human_nonCNV_mirnas = [] # store N mirnas for human non-CNV genes
    human_nonCNV_sites = []
    human_nonCNV_ratio = []
    chimp_nonCNV_hc_mirnas = [] # store N mirnas for chimp non-CNV genes with human CNV ortholog 
    chimp_nonCNV_hc_sites = []
    chimp_nonCNV_hc_ratio = []
    chimp_nonCNV_hn_mirnas = [] # store N mirnas for chimp non-CNV genes with human non-CNV ortholo
    chimp_nonCNV_hn_sites = []
    chimp_nonCNV_hn_ratio = []


    # populate the lists
    for gene in human_mirna_regulation:
        if chimp_mirna_regulation[gene][-1] == 'non-CNV':
            if human_mirna_regulation[gene][-1] == 'CNV':
                human_CNV_mirnas.append(human_mirna_regulation[gene][0])
                human_CNV_sites.append(human_mirna_regulation[gene][1])
                human_CNV_ratio.append(human_mirna_regulation[gene][2])
                chimp_nonCNV_hc_mirnas.append(chimp_mirna_regulation[gene][0])
                chimp_nonCNV_hc_sites.append(chimp_mirna_regulation[gene][1])
                chimp_nonCNV_hc_ratio.append(chimp_mirna_regulation[gene][2])
            elif human_mirna_regulation[gene][-1] == 'non-CNV':
                human_nonCNV_mirnas.append(human_mirna_regulation[gene][0])
                human_nonCNV_sites.append(human_mirna_regulation[gene][1])
                human_nonCNV_ratio.append(human_mirna_regulation[gene][2])
                chimp_nonCNV_hn_mirnas.append(chimp_mirna_regulation[gene][0])
                chimp_nonCNV_hn_sites.append(chimp_mirna_regulation[gene][1])
                chimp_nonCNV_hn_ratio.append(chimp_mirna_regulation[gene][2])
            
    return (human_CNV_mirnas, human_CNV_sites, human_CNV_ratio, 
            human_nonCNV_mirnas, human_nonCNV_sites, human_nonCNV_ratio, 
            chimp_nonCNV_hc_mirnas, chimp_nonCNV_hc_sites, chimp_nonCNV_hc_ratio, 
            chimp_nonCNV_hn_mirnas, chimp_nonCNV_hn_sites, chimp_nonCNV_hn_ratio)


def test_miRNA_regulation_CNV_non_CNV_genes(human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, paired, perry_only):
    '''
    (file, file, file, file, str, str) -> tuple
    Performs a Wilcoxon signed-rank test or a Wilcoxon sum rank test between orthologous human CNV and chimp non-CNV genes for the number of 
    miRNAs, sites and the number of sites per miRNA. Performs a Wilcoxon signed-rank test or a Wilcoxon sum rank test between orthologs non-CNV genes as control
    Return a tuple of 2-item tuple containing the z-value and the p-value
    '''

    # get the number of miRNAs, sites and sites per miRNA for human CNV genes and human and chimp non-CNV genes
    mirna_regulation = get_human_chimp_miRNA_regulation(human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only)

    human_CNV_mirnas = mirna_regulation[0]
    human_CNV_sites = mirna_regulation[1]
    human_CNV_ratio = mirna_regulation[2]
    human_nonCNV_mirnas = mirna_regulation[3]
    human_nonCNV_sites = mirna_regulation[4]
    human_nonCNV_ratio = mirna_regulation[5] 
    chimp_nonCNV_hc_mirnas = mirna_regulation[6]
    chimp_nonCNV_hc_sites = mirna_regulation[7]
    chimp_nonCNV_hc_ratio = mirna_regulation[8]
    chimp_nonCNV_hn_mirnas = mirna_regulation[9]
    chimp_nonCNV_hn_sites = mirna_regulation[10]
    chimp_nonCNV_hn_ratio = mirna_regulation[11]
    
    # compute the Wilcoxon rank sum tests
    from scipy import stats


    if paired == 'yes':
        diff_mirnas = stats.wilcoxon(human_CNV_mirnas, chimp_nonCNV_hc_mirnas)
        diff_sites = stats.wilcoxon(human_CNV_sites, chimp_nonCNV_hc_sites)
        diff_ratio = stats.wilcoxon(human_CNV_ratio, chimp_nonCNV_hc_ratio)
        diff_control_mirnas = stats.wilcoxon(human_nonCNV_mirnas, chimp_nonCNV_hn_mirnas)
        diff_control_sites = stats.wilcoxon(human_nonCNV_sites, chimp_nonCNV_hn_sites)
        diff_control_ratio = stats.wilcoxon(human_nonCNV_ratio, chimp_nonCNV_hn_ratio)

    elif paired == 'no':
        diff_mirnas = stats.ranksums(human_CNV_mirnas, chimp_nonCNV_hc_mirnas)
        diff_sites = stats.ranksums(human_CNV_sites, chimp_nonCNV_hc_sites)
        diff_ratio = stats.ranksums(human_CNV_ratio, chimp_nonCNV_hc_ratio)
        diff_control_mirnas = stats.ranksums(human_nonCNV_mirnas, chimp_nonCNV_hn_mirnas)
        diff_control_sites = stats.ranksums(human_nonCNV_sites, chimp_nonCNV_hn_sites)
        diff_control_ratio = stats.ranksums(human_nonCNV_ratio, chimp_nonCNV_hn_ratio)
        
    return (diff_mirnas, diff_sites, diff_ratio,
            diff_control_mirnas, diff_control_sites, diff_control_ratio)


def get_UTR_length_CNV_non_CNV_genes(UTR_file, human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only, species = 'human'):
    '''
    (file, file, file, file, file, str, str) -> tuple
    Returns a 4-item tuple containing the lists of UTR length CNV and non-CNV orthologous genes in human and chimp
    Option to include all human CNV genes or only human CNV genes from the same study identifying chimp CNVs
    The lists must be ordered to perform a paired test
    '''

    # get the sequence of all UTRs
    human_UTR = UTR_sequence(UTR_file, species)
    chimp_UTR = chimp_UTR_sequence(UTR_file)

    # make a set of orthologous genes
    orthologs = get_human_chimp_orthologs(human_CNV_miRNA_file, chimp_CNV_miRNA_file)

    # make a dictionnary for each gene in the CNV_miRNA file with a list containing the transcript name and the CNV status
    human_target_genes = {}
    human_cnv = open(human_CNV_miRNA_file, 'r')
    human_cnv.readline()
    for line in human_cnv:
        line = line.rstrip()
        if line != '':
            line = line.split()
            human_target_genes[line[0]] = [line[1], line[-1]]
    human_cnv.close()

    chimp_target_genes = {}
    chimp_cnv = open(chimp_CNV_miRNA_file, 'r')
    chimp_cnv.readline()
    for line in chimp_cnv:
        line = line.rstrip()
        if line != '':
            line = line.split()
            chimp_target_genes[line[0]] = [line[1], line[-1]]
    chimp_cnv.close()


    # remove genes that are not orthologs
    human_to_remove = []
    for gene in human_target_genes:
        if gene not in orthologs:
            human_to_remove.append(gene)
    for gene in human_to_remove:
        del human_target_genes[gene]

    chimp_to_remove = []
    for gene in chimp_target_genes:
        if gene not in orthologs:
            chimp_to_remove.append(gene)
    for gene in chimp_to_remove:
        del chimp_target_genes[gene]

    # use only the human CNV genes from the same study identifying chimp CNV genes or use all human CNV genes
    # if perry_only == yes : use only human and chimp CNV genes from the same study
    if perry_only == 'yes':
        perry_human_CNV_genes = get_CNV_genes_Perry_study(CNV_file, all_gene_file)
        # change the CNV status of the human target genes
        for gene in human_target_genes:
            if gene in perry_human_CNV_genes:
                human_target_genes[gene][-1] = 'CNV'
            elif gene not in perry_human_CNV_genes:
                human_target_genes[gene][-1] = 'non-CNV'

    # add the UTR length to the list of each human and chimp orthologous target gene
    for gene in human_target_genes:
        transcript = human_target_genes[gene][0]
        human_target_genes[gene].insert(-1, len(human_UTR[transcript]))
    for gene in chimp_target_genes:
        transcript = chimp_target_genes[gene][0]
        chimp_target_genes[gene].insert(-1, len(chimp_UTR[transcript]))

    # create lists to the store the UTR length of CNV and non-CNV target genes
    human_CNV_UTR = []
    human_nonCNV_UTR = []
    chimp_nonCNV_hc_UTR = []
    chimp_nonCNV_hn_UTR = []
    
    # partition the UTR length of target genes based on CNV status
    # lists must be ordered for each orthologous gene to perform a paired test
    for gene in human_target_genes:
        if chimp_target_genes[gene][-1] == 'non-CNV':
            if human_target_genes[gene][-1] == 'CNV':
                human_CNV_UTR.append(human_target_genes[gene][1])
                chimp_nonCNV_hc_UTR.append(chimp_target_genes[gene][1])
            elif human_target_genes[gene][-1] == 'non-CNV':
                human_nonCNV_UTR.append(human_target_genes[gene][1])
                chimp_nonCNV_hn_UTR.append(chimp_target_genes[gene][1])
    
    return (human_CNV_UTR, human_nonCNV_UTR, chimp_nonCNV_hc_UTR, chimp_nonCNV_hn_UTR)
       

def test_UTR_length_CNV_non_CNV_genes(UTR_file, human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only, paired, species = 'human'):
    '''
    (file, str, file, file, file, file, str) -> tuple
    Performs a Wilcoxon paired test between orthologous human CNV and chimp non-CNV genes for UTR length.
    Performs a Wilcoxon paired test between orthologs non-CNV genes as control
    Return a tuple of 2-item tuple containing the z-value and the p-value
    '''

    # get the length of UTR for CNV and non-CNV target genes
    target_UTR = get_UTR_length_CNV_non_CNV_genes(UTR_file, human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only, species = 'human')

    human_CNV_UTR = target_UTR[0]
    human_nonCNV_UTR = target_UTR[1]
    chimp_nonCNV_hc_UTR = target_UTR[2]
    chimp_nonCNV_hn_UTR = target_UTR[3]

    # performs the Wilcoxon rank sum test
    from scipy import stats
    if paired == 'yes':
        diff_UTR = stats.wilcoxon(human_CNV_UTR, chimp_nonCNV_hc_UTR)
        diff_control_UTR = stats.wilcoxon(human_nonCNV_UTR, chimp_nonCNV_hn_UTR)
    elif paired == 'no':
        diff_UTR = stats.ranksums(human_CNV_UTR, chimp_nonCNV_hc_UTR)
        diff_control_UTR = stats.ranksums(human_nonCNV_UTR, chimp_nonCNV_hn_UTR)

    return diff_UTR, diff_control_UTR


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



def print_results_human_chimp_tests(UTR_file, human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only, paired, species = 'human'):
    '''
    (file, str) -> None
    Print the results of the Wilcoxon paired tests comparing miRNA regulation and UTR length between human and chimp orthologs
    '''

    # get the list of mirnas, sites and ratio for humand and chimp orthologs
    mirna_regulation = get_human_chimp_miRNA_regulation(human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only)

    # performs paired tests between humand and chimp orthologs for mirna regulation
    test_regulation = test_miRNA_regulation_CNV_non_CNV_genes(human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, paired, perry_only)

    # get the list of UTR length for humand and chimp orthologs
    mirna_UTR = get_UTR_length_CNV_non_CNV_genes(UTR_file, human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only, species = 'human')

    # performs paired tests between human and chimp orthologs for UTR length
    test_UTR = test_UTR_length_CNV_non_CNV_genes(UTR_file, human_CNV_miRNA_file, chimp_CNV_miRNA_file, CNV_file, all_gene_file, perry_only, paired, species = 'human')

   
    wilcoxon_regulation = ['mirnas_human_cnv_chimp_noncnv:', 'sites_human_cnv_chimp_noncnv:', 'ratio_human_cnv_chimp_noncnv:',
                           'mirnas_human_noncnv_chimp_noncnv:', 'sites_human_noncnv_chimp_noncnv:', 'ratio_human_noncnv_chimp_noncnv:']

    if paired == 'yes':
        print('Wilcoxon signed rank tests for differences between humand and chimp in miRNA regulation:')
    elif paired == 'no':
        print('Wilcoxon sum rank tests for differences between humand and chimp in miRNA regulation:')
                
    for i in range(len(test_regulation)):
        print(wilcoxon_regulation[i], 'z-score = ', test_regulation[i][0], 'p = ', test_regulation[i][1], sep = '\t')
            
    print('\n')
    print('N human CNV genes:' + '\t' + str(len(mirna_regulation[0])))
    print('N human non-CNV genes:' + '\t' + str(len(mirna_regulation[3])))
    print('\n')
    print('\t' + 'mean' + '\t' + 'std_error')
    headers = ['human_CNV_mirnas', 'human_CNV_sites', 'human_CNV_ratio',
               'human_nonCNV_mirnas', 'human_nonCNV_sites', 'human_nonCNV_ratio',
               'chimp_nonCNV_human_cnv_mirnas', 'chimp_nonCNV_human_cnv_sites', 'chimp_nonCNV_human_cnv_ratio',
               'chimp_nonCNV_human_noncnv_mirnas', 'chimp_nonCNV_human_noncnv_sites', 'chimp_nonCNV_human_noncnv_ratio']
      
    for i in range(len(mirna_regulation)):
        mean_stderr = compute_mean_std_error(mirna_regulation[i])
        print(headers[i], mean_stderr[0], mean_stderr[1], sep = '\t')

    print('\n')
        
    if paired == 'yes':
        print('Wilcoxon signed rank tests for differences between humand and chimp in UTR length:')
    elif paired == 'no':
        print('Wilcoxon sum rank tests for differences between humand and chimp in UTR length:')
   
    print('human CNV - chimp nonCNV:', 'z-score:', test_UTR[0][0], 'p', test_UTR[0][1], sep = '\t')
    print('human nonCNV - chimp nonCNV:', 'z-score:', test_UTR[1][0], 'p', test_UTR[1][1], sep = '\t')

    print('\n')

    print('\t' + 'mean' + '\t' + 'std_error')
    utr_header = ['human_CNV_UTR', 'human_nonCNV_UTR', 'chimp_nonCNV_human_cnv_UTR', 'chimp_nonCNV_human_noncnv_UTR']

    for i in range(len(mirna_UTR)):
        mean_stderr = compute_mean_std_error(mirna_UTR[i])
        print(utr_header[i], mean_stderr[0], mean_stderr[1], sep = '\t')


