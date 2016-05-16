def transcript_gene_pairs(summary_counts):
    '''
    (file) -> dict
    Return a dictionnary with the transcripts as key and the corresponding gene as value
    '''
    # open the file, create a dictionnary 
    summary = open(summary_counts, 'r')
    header = summary.readline()
    transcripts = {}

    # populate the dictionnary with the transcript: gene pairs
    for line in summary:
        line = line.rstrip()
        if line != '':
            line = line.split()
            transcripts[line[0]] = line[1]

    summary.close()
    return transcripts



def fly_transcript_gene_pairs(UTR_file):
    '''
    (file) -> dict
    Return a dictionnary with the transcripts as key and the corresponding gene as value
    Precondition: species = drosophila
    '''

    # in drosophila, the summary counts includes only transcripts and not genes
    # there is a single transcript per gene in droso
    utr = open(UTR_file, 'r')
    header = utr.readline()
    transcripts = {}

    # populate the dictionnary with the transcript: gene pairs
    for line in utr:
        if '\t7227\t' in line:
            line = line.rstrip().split()
            transcripts[line[2]] = line[1]
    utr.close()
    return transcripts
    

def UTR_sequence(UTR_file, species):
    '''
    (file, str) -> dict
    Return a dictionnary with the UTR sequence of each transcript for a given species
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
            if species == 'human':
                if '9606' in line:
                    utr_sequences[transcript] = sequence
            elif species == 'worm':
                if '6239' in line:
                    utr_sequences[transcript] = sequence
            elif species == 'zebrafish':
                if '7955' in line:
                    utr_sequences[transcript] = sequence
            elif species == 'fly':
                if '7227' in line:
                    transcript = line[2]
                    utr_sequences[transcript] = sequence

    utr.close()
    return utr_sequences


def miRNA_family(miRNA_family_info, species):
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
            if species == 'human':
                if line[2] == '9606':
                    if line[1] in family:
                        family[line[1]].append(line[3])
                    else:
                        family[line[1]] = [line[3]]
            elif species == 'worm':
                if line[2] == '6239':
                    if line[1] in family:
                        family[line[1]].append(line[3])
                    else:
                        family[line[1]] = [line[3]]
            elif species == 'zebrafish':
                if line[2] == '7955':
                    if line[1] in family:
                        family[line[1]].append(line[3])
                    else:
                        family[line[1]] = [line[3]]
            elif species == 'fly':
                if line[2].startswith('dme'):
                    if line[1] in family:
                        family[line[1]].append(line[2])
                    else:
                        family[line[1]] = [line[2]]

    miRNA.close()
    return family


def human_CNV_genes(CNV_file, all_gene_file):
    '''
    (file, file) -> set
    Returns a set of CNV genes
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
            if line[4] == 'CNV':
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


def human_CNV_genes_single_study(CNV_file, all_gene_file, study):
    '''
    (file, file, study) -> set
    Returns a set of CNV genes for a given study in the CNV_file, eliminating genes not found in the all_gene_file
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
            if line[4] == 'CNV' and line[6] == study:
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


def worm_CNV_genes(CNV_file, all_gene_file):
    '''
    (file, file) -> set
    Return a set of CNV genes
    '''

    # get a set of valid worm genes
    valid_genes = sort_valid_worm_genes(all_gene_file)

    # make a set of CNV genes
    CNV = open(CNV_file, 'r')
    CNV_genes = set()
    header = CNV.readline()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            line = line.split()
            CNV_genes.add(line[1])

    # remove non valid CNV genes
    to_remove = []
    for gene in CNV_genes:
        if gene not in valid_genes:
            to_remove.append(gene)
    for gene in to_remove:
        CNV_genes.discard(gene)

    CNV.close()
    return CNV_genes



def zebrafish_transcript_gene_pairs(ensembl_gene_transcripts):
    '''
    (file) -> dict
    Return a dictionnary with ensembl transcript ID as key and ensembl gene ID as value
    '''
    
    # make a dictionnary with the ensembl transcript : gene pairs
    transcripts_genes = {}
    transcripts = open(ensembl_gene_transcripts, 'r')
    transcripts.readline()
    for line in transcripts:
        line = line.rstrip()
        if line !='':
            line = line.split()
            transcripts_genes[line[1]] = line[0]

    transcripts.close()
    return transcripts_genes


def zebrafish_valid_chromosomes():
    '''
    () -> set
    Return a set of zebrafish chromosome IDs from chr1 to chr25
    '''

    # make a set of valid chromosomes
    valid_chromosomes = set()
    for i in range(1, 26):
        valid_chromosomes.add('chr' + str(i))

    return valid_chromosomes

            
def zebrafish_CNV_genes(CNV_positions_file, ensembl_transcript_coordinates, ensembl_gene_transcripts):
    '''
    (file, file, file) -> set
    Return a set of zebrafish CNV genes
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
            CNV_positions[i] = ['chr' + line[0], int(line[1]), int(line[2])]
            i += 1

    # make a dictionnary with the ensembl transcript : gene pairs
    transcripts_genes = zebrafish_transcript_gene_pairs(ensembl_gene_transcripts)
    
    # make a set of valid chromosomes
    valid_chromosomes = zebrafish_valid_chromosomes()
        
    # make a dictionnary to store the coordinates of the ensembl transcripts {TS1 : [chromo, start, end]}
    transcripts_coord = open(ensembl_transcript_coordinates, 'r')
    transcripts_coord.readline()
    transcripts_positions = {}
    for line in transcripts_coord:
        line = line.rstrip()
        if line !='':
            line = line.split()
            if line[2] in valid_chromosomes:
                transcripts_positions[line[1]] = [line[2], int(line[4]), int(line[5])]

    print(len(transcripts_positions))
    
    # search for CNV transcripts using the transcript and CNV coordinates and store in a set
    # report transcripts that overlap with a CNV even in the UTR or non-coding regions
    done = 0
    zebrafish_CNV_genes = set()
    for transcript in transcripts_positions:
        if done % 5 == 0:
            print(done, len(zebrafish_CNV_genes), sep = '\t')
        done += 1
        if transcripts_genes[transcript] not in zebrafish_CNV_genes:     # no need to search again if gene already CNV
            for CNV in CNV_positions:
                if transcripts_positions[transcript][0] == CNV_positions[CNV][0]:
                    ts_pos = set(range(transcripts_positions[transcript][1], transcripts_positions[transcript][2] + 1))
                    cnv_pos = set(range(CNV_positions[CNV][1], CNV_positions[CNV][2] + 1))
                    # if CNV and transcript coordinates overlap keep transcript
                    if len(ts_pos.intersection(cnv_pos)) != 0:
                        zebrafish_CNV_genes.add(transcripts_genes[transcript]) 
                        break                                            # no need to search again if gene already CNV
        
    cnv_coord.close()
    transcripts_coord.close()
    return zebrafish_CNV_genes



def fly_CNV_genes(CNV_positions_file, fly_gene_coordinates):
    '''
    (file, file) -> dict
    Returns a dictionnary of fly CNV genes with the FlyBase gene_ID as key and gene name / symbol as value
    '''

    # make a dictionnary to store the coordinates of each CNV {CNV1:[chromo, start, end]}
    cnv_coord = open(CNV_positions_file, 'r')
    cnv_coord.readline()
    CNV_positions = {}
    i = 0
    for line in cnv_coord:
        line = line.rstrip()
        if line != '':
            line = line.split()
            CNV_positions[i] = [line[2], int(line[0]), int(line[1])]
            i += 1

    # make a dictionnary with the flybase ID : gene name pairs
    gene_IDs = {}

    # make a dictionnary to store the gene coordinates using the flybase gene_ID as key {FBgn1 : [chromo, start, end]}
    gene_positions = {}
    
    # make a set of valid chromosomes
    valid_chromosomes = {'2R', '3R', '4', '3L', '2L', 'X'}
    
    # open and read file
    all_genes = open(fly_gene_coordinates, 'r')
    for line in all_genes:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'FlyBase' and line[2] == 'gene' and line[0] in valid_chromosomes:
                gene = line[8][line[8].index('ID')+3:line[8].index(';')]
                name = line[8][line[8].index('Name')+5: line[8].index(';', line[8].index('Name')+5)]
                gene_IDs[gene] = name
                start = int(line[3])
                end = int(line[4])
                gene_positions[gene] = [line[0], start, end]
         
    print(len(gene_positions))
        
    # search for CNV genes using the gene and CNV coordinates and store in a set
    # report transcripts that overlap with a CNV even in the UTR or non-coding regions
    done = 0
    fly_CNV_genes = set()
    for gene in gene_positions:
        if done % 5 == 0:
            print(done, len(fly_CNV_genes), sep = '\t')
        done += 1
        for CNV in CNV_positions:
            if gene_positions[gene][0] == CNV_positions[CNV][0]:
                gene_pos = set(range(gene_positions[gene][1], gene_positions[gene][2] + 1))
                cnv_pos = set(range(CNV_positions[CNV][1], CNV_positions[CNV][2] + 1))
                # if CNV and gene coordinates overlap keep gene
                if len(gene_pos.intersection(cnv_pos)) != 0:
                    fly_CNV_genes.add(gene)
                    break   # no need to search again if gene already CNV


    # remove genes that are not CNV
    to_remove = []
    for gene in gene_IDs:
        if gene not in fly_CNV_genes:
            to_remove.append(gene)
    for gene in to_remove:
        del gene_IDs[gene]

    cnv_coord.close()
    all_genes.close()
    return gene_IDs



def sort_valid_worm_genes(all_gene_file):
    '''
    (file) -> set
    Return a set of valid C. elegans genes that are alive in the current genome annotation file
    '''

    # make a set containing gene id, symbol and gene names
    valid_genes = set()
    IDs = open(all_gene_file, 'r')
    for line in IDs:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            if line[-1] == 'Live':
                while '' in line:
                    line.remove('')
                for item in line:
                    valid_genes.add(item)

    IDs.close()
    return valid_genes


def sort_valid_human_genes(all_gene_file):
    '''
    (file) -> set
    Return a set of valid human genes
    '''

    human_genes = set()
    all_genes = open(all_gene_file, 'r')
    header = all_genes.readline()
    for line in all_genes:
        line = line.rstrip()
        if line != '':
            line = line.split()
            human_genes.add(line[2])

    all_genes.close()
    return human_genes


def sort_valid_zebrafish_genes(ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only = 'yes'):
    '''
    (file, file) -> set
    Return a set of valid zebrafish ensembl gene names
    By default: do not consider transcripts not assigned to chromosomes because CNVs are reported only for chromosomes
    '''

    if chromosome_only == 'yes':
        # make a dictionnary of transcript : genes pairs
        transcripts = zebrafish_transcript_gene_pairs(ensembl_gene_transcripts)
        
        # make a set of valid chromosomes
        valid_chromosome = zebrafish_valid_chromosomes()
        
        # make a set of transcripts located on chromosomes
        valid_transcripts = set()
        transcripts_coord = open(ensembl_transcript_coordinates, 'r')
        transcripts_coord.readline()
        for line in transcripts_coord:
            line = line.rstrip()
            if line != '':
                line = line.split()
                if line[2] in valid_chromosome:
                    valid_transcripts.add(line[1])

        # remove transcripts from the dictionnary that are not on chromosomes
        to_remove = []
        for transcript_ID in transcripts:
            if transcript_ID not in valid_transcripts:
                to_remove.append(transcript_ID)
        for transcript_ID in to_remove:
            del transcripts[transcript_ID]

        # make a set of gene names
        zebrafish_genes = set()
        for transcript_ID in transcripts:
            zebrafish_genes.add(transcripts[transcript_ID])

    elif chromosome_only == 'no':
        # make a set of gene names including genes in scaffold
        zebrafish_genes = set()
        all_genes = open(ensembl_gene_transcripts, 'r')
        header = all_genes.readline()
        for line in all_genes:
            line = line.rstrip()
            if line != '':
                line = line.split()
                zebrafish_genes.add(line[0])

    if chromosome_only == 'yes':
        transcripts_coord.close()
    elif chromosome_only == 'no':
        all_genes.close()
    return zebrafish_genes




def sort_valid_fly_genes(fly_gene_coordinates):
    '''
    (file) -> set
    Returns a set of FlyBase genes from the gff file fly_gene_coordinates, keeping both the FlyBase gene ID and the gene name in the set
    '''

    # make a set of valid chromosomes:
    valid_chromosomes = {'2R', '3R', '4', '3L', '2L', 'X'}

    # make a set of valid_genes
    valid_genes = set()
    
    # open and read file
    all_genes = open(fly_gene_coordinates, 'r')

    for line in all_genes:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'FlyBase' and line[2] == 'gene' and line[0] in valid_chromosomes:
                gene = line[8][line[8].index('ID')+3:line[8].index(';')]
                name = line[8][line[8].index('Name')+5: line[8].index(';', line[8].index('Name')+5)]
                valid_genes.add(gene)
                valid_genes.add(name)

    all_genes.close()
    return valid_genes


def human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites):
    '''
    (file, str, str) -> dict
    Return a dictionnary with the number of miRNA regulators, miRNA binding sites and number of sites per miRNA for each target transcript
    of a given species exctracted from the summary_count file. Use either only conserved sites or all sites.
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
            if species == 'human': 
                if line[3] == '9606':
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
    

def worm_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info):
    '''
    (file, str, str, file) -> dict
    Return a dictionnary with the number of miRNA regulators, miRNA binding sites and number of sites per miRNA for each target transcript
    exctracted from the summary_count file. Use either only conserved sites or all sites, but use only miRNAs in the miRNA_family_info file
    (ie. does not take into account miR* present in the summary_count file)
    '''

    # open the summary_count file for reading
    summary = open(summary_counts, 'r')
    header = summary.readline()

    # create a dictionnary to store the miRNA family and number of sites for each transcripts
    # {transcript_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2]]}
    regulated_transcripts = {}

    # create a dictionary with the seeds as keys and list of miRNAs as value
    seeds = miRNA_family(miRNA_family_info, species)

    # read the summary_count file
    for line in summary:
        line = line.rstrip()
        if line != '':
            line = line.split()
            transcript = line[0]
            family = line[2]
            N_conserved_sites = int(line[4])
            N_poorly_conserved_sites = int(line[8])
            if species == 'worm':
                # do not take into account targets of star miRNAs
                if line[3] == '6239' and family in seeds:
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


def zebrafish_miRNA_regulation_TargetScan_sites(summary_counts, species, miRNA_family_info):
    '''
    (file, str, file) -> dict
    Return a dictionnary with the number of miRNA regulators, miRNA binding sites and number of sites per miRNA for each target transcript
    exctracted from the summary_count file
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
            if species == 'zebrafish': 
                if '7955' in line:
                    family = line[line.index('7955') -1] # the index of N_sites is not always the same but is always before species_ID
                    N_sites = int(line[line.index('7955') + 1]) # the index of N_sites is not always the same but is always after species_ID
                    if transcript in regulated_transcripts:
                        regulated_transcripts[transcript].append([family, N_sites])
                    else:
                        regulated_transcripts[transcript] = [[family, N_sites]]
                                    
    summary.close()
    return regulated_transcripts


def fly_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info):
    '''
    (file, str, str, file) -> dict
    Return a dictionnary with the number of miRNA regulators, miRNA binding sites and number of sites per miRNA for each target transcript
    exctracted from the summary_count file. Use either only conserved sites or all sites
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
            transcript = line[2]
            family = line[0]
            N_conserved_sites = int(line[3])
            N_poorly_conserved_sites = int(line[7])
            if species == 'fly':
                if line[1] == '7227':
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


def human_worm_make_miRNA_regulators_table(summary_counts, species, conservation_sites, miRNA_family_info, CNV_file, all_gene_file, outputfile):
    '''
    (file, str, str, file, file, file, file) -> file
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not
    '''

    # get the gene ID for each target transcript
    transcripts = transcript_gene_pairs(summary_counts)

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

    # get the CNV status and sort the valid genes
    # get the miRNA regulators and binding sites
    if species == 'human':
        CNV_genes = human_CNV_genes(CNV_file, all_gene_file)
        valid_genes = sort_valid_human_genes(all_gene_file)
        regulators = human_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites)
    elif species == 'worm':
        CNV_genes = worm_CNV_genes(CNV_file, all_gene_file)
        valid_genes = sort_valid_worm_genes(all_gene_file)
        regulators = worm_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)
        
    # create a dictionnary to store info about miRNA regultors for each gene
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, CNV_status]}
    gene_regulated = {}
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


def zebrafish_make_miRNA_regulators_table(summary_counts, species, miRNA_family_info, CNV_genes_file, ensembl_gene_transcripts, ensembl_transcript_coordinates, outputfile):
    '''
    (file, str, file, file, file, file, file) -> file
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites,
    the mean number of sites / mirna and whether the gene is in a CNV or not
    '''

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

    # grab the CNV genes
    CNV = open(CNV_genes_file, 'r')
    CNV_genes = set()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            CNV_genes.add(line)

    # sort the valid genes
    valid_genes = sort_valid_zebrafish_genes(ensembl_gene_transcripts, ensembl_transcript_coordinates, chromosome_only = 'yes')
    
    # get the miRNA regulators and binding sites
    regulators = zebrafish_miRNA_regulation_TargetScan_sites(summary_counts, species, miRNA_family_info)

    # create a dictionnary to store info about miRNA regultors for each gene
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, CNV_status]}
    gene_regulated = {}
    for transcript_ID in regulators:
        gene = transcript_ID[:transcript_ID.index('.')] # the transcript ID is not a valid ensembl ID but the gene_ID + .number
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

    CNV.close()
    newfile.close()


def fly_make_miRNA_regulators_table(UTR_file, miRNA_family_info, CNV_genes_file, species, fly_gene_coordinates, summary_counts, conservation_sites, outputfile):
    '''
    (file, file, file, str, file, file, str, file) -> file
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites (conserved or all),
    the mean number of sites / mirna and whether the gene is in a CNV or not
    '''

    # get the gene ID for each target transcript
    transcripts = fly_transcript_gene_pairs(UTR_file)

    # get the list of miRNAs for each family
    family = miRNA_family(miRNA_family_info, species)

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

    # get the valid genes
    valid_genes = sort_valid_fly_genes(fly_gene_coordinates)

    # get the mirna info for each transcript
    regulators = fly_miRNA_regulation_TargetScan_sites(summary_counts, species, conservation_sites, miRNA_family_info)

    # create a dictionnary to store info about miRNA regultors for each gene
    # using a single transcript per gene
    # {gene_1: [[family_1, N_sites_family_1], [family_2, N_sites_family_2], transcript_ID, CNV_status]}
    gene_regulated = {}
    for transcript_ID in regulators:
        gene = transcripts[transcript_ID]
        if gene not in gene_regulated:
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
        newfile.write(gene + '\t')
        newfile.write(gene_regulated[gene][-2] + '\t')
        N_mirnas = 0
        N_sites = 0
        for pair in gene_regulated[gene][:-2]:
            N_mirnas += len(family[pair[0]])
            N_sites += pair[1]
        newfile.write(str(N_mirnas) + '\t' + str(N_sites) + '\t' + str(N_sites/ N_mirnas) + '\t' + gene_regulated[gene][-1] + '\n')

    newfile.close()

