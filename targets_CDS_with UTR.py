from CNV_miRNAs_TargetScan import *
from targetscan_same_size_UTR import *
from scipy import stats



def merge_human_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, all_gene_file, species):
    '''
    (file, file, file, file, str) -> dict
    Return a dictionnary for genes in the CNV_mirna file with merged CDS
    from the CDS file and the 3' UTR file from the UTR file
    '''

    # get the 3'UTR sequences for all transcripts
    transcripts_UTR = UTR_sequence(UTR_file, species)

    # make a dictionnary with the gene : transcript pairs from the CNV mirna file
    genes_to_transcripts = {}

    # open CNV miRNA file for reading
    cnvs = open(CNV_miRNA_file, 'r')
    # skipe header
    cnvs.readline()
    # populate dict
    for line in cnvs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            genes_to_transcripts[line[0]] = line[1]
    # close line after reading
    cnvs.close()

    # make a dictionnary with the ensembl gene : CDS sequence from the CDS file
    genes_CDS = {}

    # open CDS file for reading
    cds = open(CDS_file, 'r')
    for line in cds:
        line = line.rstrip()
        if line.startswith('>'):
            line = line[1:].split(' ')
            name = line[3][line[3].index('gene:') + 5:]
            genes_CDS[name] = ''
        else:
            genes_CDS[name] += line
    # close file after reading
    cds.close()
            
    # replace Ts with Us
    for gene in genes_CDS:
        genes_CDS[gene] = genes_CDS[gene].replace('T', 'U')

    # make a dictionnary with the gene: UTR pairs
    genes_UTR = {}

    # populate dict with the gene : 3' UTR pairs
    for gene in genes_to_transcripts:
        genes_UTR[gene] = transcripts_UTR[genes_to_transcripts[gene]]

    # make a dictionnary to hold the gene : Ensembl gene name pairs
    genes_to_ensembl = {}

    # open the all gene file for reading
    all_genes = open(all_gene_file, 'r')
    all_genes.readline()

    # populate dictionnary with gene : ensembl name pairs
    for line in all_genes:
        line = line.rstrip()
        if line != '':
            line = line.split()
            for i in range(len(line)):
                if 'Ensembl' in line[i]:
                    gene = line[2]
                    ensembl_name = line[i][line[i].index('Ensembl:') + 8 : line[i].index('Ensembl:') + 8+ 15]
                    genes_to_ensembl[gene] = ensembl_name
    # close file after reading
    all_genes.close()

    # make a dictionary with the CDS and UTR sequences merged
    genes_CDS_UTR = {}
     
    # go though the gene : UTR dict and apopulate the dict with the CDS and the UTR
    for gene in genes_UTR:
        if gene in genes_to_ensembl:
            ensembl = genes_to_ensembl[gene]
            if ensembl in genes_CDS:
                CDS_seq = genes_CDS[ensembl]
                genes_CDS_UTR[gene] = CDS_seq + genes_UTR[gene]

    return genes_CDS_UTR

    
def prepare_targetscan_UTR_input(CNV_miRNA_file, UTR_file, CDS_file, all_gene_file, species, outputfile):
    '''
    (file, file, file, file, str, file) -> file
    Save to file the CDS + UTR sequence as input for targetscan site prediction
    '''

    # merged the CDS and UTR sequence
    if species == 'human':
        genes_CDS_UTR = merge_human_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, all_gene_file, species)
    elif species == 'worm':
        genes_CDS_UTR = merge_worm_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, all_gene_file, species)
    elif species == 'fly':
        genes_CDS_UTR = merge_fly_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, species)
    elif species == 'zebrafish':
        genes_CDS_UTR = merge_zebrafish_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, all_gene_file, species)
        
    # open file for writing
    newfile = open(outputfile, 'w')

    # write content of dict to file
    for gene in genes_CDS_UTR:
        if species == 'human':
            newfile.write(gene + '\t' + '9606' + '\t' + genes_CDS_UTR[gene] + '\n')
        elif species == 'worm':
            newfile.write(gene + '\t' + '6239' + '\t' + genes_CDS_UTR[gene] + '\n')
        elif species == 'fly':
            newfile.write(gene + '\t' + '7227' + '\t' + genes_CDS_UTR[gene] + '\n')
        elif species == 'zebrafish':
            newfile.write(gene + '\t' + '7955' + '\t' + genes_CDS_UTR[gene] + '\n')
                        
    # close file after writing
    newfile.close()
    
    
def merge_worm_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, all_gene_file, species):
    '''
    (file, file, file, file, str) -> dict
    Return a dictionnary for genes in the CNV_mirna file with merged CDS
    from the CDS file and the 3' UTR file from the UTR file
    '''

    # get the 3'UTR sequences for all transcripts
    transcripts_UTR = UTR_sequence(UTR_file, species)

    # make a dictionnary with the gene : transcript pairs from the CNV mirna file
    genes_to_transcripts = {}

    # open CNV miRNA file for reading
    cnvs = open(CNV_miRNA_file, 'r')
    # skipe header
    cnvs.readline()
    # populate dict
    for line in cnvs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            genes_to_transcripts[line[0]] = line[1]
    # close line after reading
    cnvs.close()

    # make a dictionnary with the WBGgene : CDS sequence from the CDS file
    genes_CDS = {}

    # open CDS file for reading
    cds = open(CDS_file, 'r')
    for line in cds:
        line = line.rstrip()
        if line.startswith('>') and 'gene' in line:
            name = line[line.index('gene=') + 5:]
            genes_CDS[name] = ''
        else:
            genes_CDS[name] += line
    # close file after reading
    cds.close()
            
    # replace Ts with Us
    for gene in genes_CDS:
        genes_CDS[gene] = genes_CDS[gene].upper().replace('T', 'U')

    # make a dictionnary with the gene: UTR pairs
    genes_UTR = {}

    # populate dict with the gene : 3' UTR pairs
    for gene in genes_to_transcripts:
        genes_UTR[gene] = transcripts_UTR[genes_to_transcripts[gene]]


    # make a dictionnary to hold the gene WBGene : gene names pairs
    WBGgene_to_genes = {}

    # open the all gene file for reading
    all_genes = open(all_gene_file, 'r')
    
    for line in all_genes:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            # keep only live genes
            if line[-1] == 'Live':
                # remove empty entries
                while '' in line:
                    line.remove('')
                # remove the gene status
                line.remove('Live')
                # remove the species ID
                line.remove('6239')
                # use the WBGgene as key and a list of genes
                WBGgene_to_genes[line[0]] = []
                for i in range(1, len(line)):
                    WBGgene_to_genes[line[0]].append(line[i])
    # close file after reading
    all_genes.close()


    # reverse the WBGgene_to_genes to make a gene_to_WBGgene dict
    gene_to_WBGgene = {}
    for WBGgene in WBGgene_to_genes:
        if len(WBGgene_to_genes[WBGgene]) == 1:
            gene = WBGgene_to_genes[WBGgene][0]
            gene_to_WBGgene[gene] = WBGgene
        elif len(WBGgene_to_genes[WBGgene]) > 1:
            for gene in WBGgene_to_genes[WBGgene]:
                gene_to_WBGgene[gene] = WBGgene


    # make a dictionary with the CDS and UTR sequences merged
    genes_CDS_UTR = {}

    # make a set of already added WBGgenes to make sure that genes with multiple names are not added multiple times
    
     
    # go though the gene : UTR dict and populate the dict with the CDS and the UTR
    for gene in genes_UTR:
        if gene in gene_to_WBGgene:
            WBGgene = gene_to_WBGgene[gene]
            if WBGgene in genes_CDS:
                CDS_seq = genes_CDS[WBGgene]
                genes_CDS_UTR[gene] = CDS_seq + genes_UTR[gene].upper().replace('T', 'U')


    return genes_CDS_UTR


def merge_fly_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, species):
    '''
    (file, file, file, file, str) -> dict
    Return a dictionnary for genes in the CNV_mirna file with merged CDS
    from the CDS file and the 3' UTR file from the UTR file
    '''

    # get the 3'UTR sequences for all transcripts
    transcripts_UTR = UTR_sequence(UTR_file, species)

    # make a set of target CNV transcripts from the CNV_miRNA_file
    target_transcripts = set()

    # open CNV miRNA file for reading
    cnvs = open(CNV_miRNA_file, 'r')
    # skipe header
    cnvs.readline()
    # populate dict
    for line in cnvs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            target_transcripts.add(line[1])
    # close line after reading
    cnvs.close()

    # make a dictionnary with the transcript : CDS sequence from the CDS file
    genes_CDS = {}

    # open CDS file for reading
    cds = open(CDS_file, 'r')
    for line in cds:
        line = line.rstrip()
        if line.startswith('>'):
            name = line[1:line.index(' ')]
            genes_CDS[name] = ''
        else:
            genes_CDS[name] += line
    # close file after reading
    cds.close()
            
    # replace Ts with Us
    for gene in genes_CDS:
        genes_CDS[gene] = genes_CDS[gene].upper().replace('T', 'U')

    # make a dictionary with the CDS and UTR sequences merged
    genes_CDS_UTR = {}

    # go though the gene : UTR dict and populate the dict with the CDS and the UTR
    for gene in target_transcripts:
        if gene in genes_CDS:
            CDS_seq = genes_CDS[gene]
            genes_CDS_UTR[gene] = CDS_seq + transcripts_UTR[gene].upper().replace('T', 'U')

    return genes_CDS_UTR



def merge_zebrafish_cds_UTR(CNV_miRNA_file, UTR_file, CDS_file, all_gene_file, species):
    '''
    (file, file, file, file, str) -> dict
    Return a dictionnary for genes in the CNV_mirna file with merged CDS
    from the CDS file and the 3' UTR file from the UTR file
    '''

    # get the 3'UTR sequences for all transcripts
    transcripts_UTR = UTR_sequence(UTR_file, species)

    # make a dictionnary to hold the transcript : gene pairs
    transcripts_to_genes = {}
    
    # open CNV miRNA file for reading
    cnvs = open(CNV_miRNA_file, 'r')
    # skipe header
    cnvs.readline()
    # populate dict
    for line in cnvs:
        line = line.rstrip()
        if line != '':
            line = line.split()
            transcripts_to_genes[line[1]] = line[0]
    # close line after reading
    cnvs.close()

    # make a dictionnary with the gene : UTR pairs
    genes_UTR = {}
    for transcript in transcripts_UTR:
        if transcript in transcripts_to_genes:
            genes_UTR[transcripts_to_genes[transcript]] = transcripts_UTR[transcript]

    # make a dictionnary to hold the gene : ensembl_transcript pairs
    genes_to_ensembl = zebrafish_transcript_gene_pairs(all_gene_file)

    # reverse the dictionnary
    ensembl_to_genes = {}
    for gene in genes_to_ensembl:
        ensembl_to_genes[genes_to_ensembl[gene]] = gene

    # make a dictionnary with the transcript : CDS sequence from the CDS file
    genes_CDS = {}

    # open CDS file for reading
    cds = open(CDS_file, 'r')
    for line in cds:
        line = line.rstrip()
        if line.startswith('>'):
            name = line[1:line.index(' ')]
            genes_CDS[name] = ''
        else:
            genes_CDS[name] += line
    # close file after reading
    cds.close()
            
    # make a dictionary with the CDS and UTR sequences merged
    genes_CDS_UTR = {}

    # go though the gene : UTR dict and populate the dict with the CDS and the UTR
    for gene in genes_UTR:
        if gene in ensembl_to_genes:
            ensembl_gene = ensembl_to_genes[gene]
            if ensembl_gene in genes_CDS:
                CDS_seq = genes_CDS[ensembl_gene]
                genes_CDS_UTR[gene] = CDS_seq.upper().replace('T', 'U') + genes_UTR[gene].upper().replace('T', 'U')

    return genes_CDS_UTR




def make_CDS_UTR_table_CNV_miRNA(targetscan_output, CNV_miRNA_file, miRNAfamily_info, species, outputfile):
    '''
    (file, file) -> file
    Extract the number of miRNAs and the number of sites from the targetscan_output file and save in outputfile
    '''

    # grab all the mirnas for a same family
    same_family = mirnas_per_family(miRNAfamily_info, species)

    # grab the seed : family name pair
    family = mirna_family_seed_name(miRNAfamily_info, species)
    # inverse the dictionnary to make a dcit with the family name as key and the seed as value
    family_name = {}
    for seed in family:
        name = family[seed][0]
        family_name[name] = seed


    # open miRNAs CNV file for reading
    CNVs = open(CNV_miRNA_file, 'r')
    # skip header
    CNVs.readline()



    if species == 'fly':
        # make a dictionnary with the transcript as key and a list with the gene name and CNV status as value
        CNV_genes = {}
        # go through the file and populate the dict
        for line in CNVs:
            line = line.rstrip()
            if line != '':
                line  = line.split()
                gene = line[0]
                transcript = line[1]
                status = line[-1]
                CNV_genes[transcript] = [gene, status]
    else:
        # make a dictionnary with the gene name as key and a list with the transcript name and CNV status as value
        CNV_genes = {}
        # go through the file and populate the dict
        for line in CNVs:
            line = line.rstrip()
            if line != '':
                line  = line.split()
                gene = line[0]
                transcript = line[1]
                status = line[-1]
                CNV_genes[gene] = [transcript, status]
    # close file after reading
    CNVs. close()

    # make a dictionnary with gene as key and a list with the sets of mirna families amd set of sites as value
    # site is defined by 'start:end:site_type'
    target_sites = {}
    # open targetscan output for reading
    targetscan = open(targetscan_output, 'r')
    targetscan.readline()
    for line in targetscan:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # determine the site
            gene = line[0]
            mirfam = line[1]
            start = line[3]
            end = line[4]
            site_type = line[8]
            site = ':'.join([start, end, site_type])
            # check if gene in dict, and add the mirfam and site as value
            if gene in target_sites:
                target_sites[gene][0].add(mirfam)
                target_sites[gene][1].add(site)
            else:
                target_sites[gene] = [{mirfam,}, {site,}]
    # close file after reading
    targetscan.close()

    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
    # write content to file
    for gene in target_sites:
        mirna_count = 0
        for mirfam in target_sites[gene][0]:
            seed = family_name[mirfam]
            mirna_count += len(same_family[seed])
        site_count = len(target_sites[gene][1])
        if species == 'fly':
            newfile.write(CNV_genes[gene][0] + '\t' + gene + '\t' + str(mirna_count) + '\t' + str(site_count) + '\t')
        else:
            newfile.write(gene + '\t' + CNV_genes[gene][0] + '\t' + str(mirna_count) + '\t' + str(site_count) + '\t')
        if mirna_count != 0:
            newfile.write(str(site_count / mirna_count) + '\t')
        else:
            newfile.write('NA' + '\t')
        newfile.write(CNV_genes[gene][1] + '\n')
    # close file ater writing
    newfile.close()





############################################



def make_CDSUTR_length_miRNA_regulation_lists(targetscan_CDS_UTR_input_file, species, CNV_miRNA_file):
    '''
    (file, str, file) -> list
    Returns a list of lists containing the length of CDS+UTR targets from the targetscan input file
    and in the same order the number of miRNAs, binding sites and ratio extracted from the CNV_miRNA_file
    '''

    # get the sequence of all CDS + UTRs for the given species
    CDS_UTR = {}

    # open file for reading
    cds = open(targetscan_CDS_UTR_input_file, 'r')
    for line in cds:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            seq = line[-1]
            CDS_UTR[gene] = len(seq)
    #close file after reading
    cds.close()

    # make a dictionnary with information for each target transcript
    targets = {}
    CNV_miRNA = open(CNV_miRNA_file, 'r')
    CNV_miRNA.readline()
    for line in CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if species == 'fly':
                # use the transcript as key
                targets[line[1]] = [int(line[2]), int(line[3]), float(line[4])]
            else:
                # use the gene as key
                targets[line[0]] = [int(line[2]), int(line[3]), float(line[4])]
    # close after reading
    CNV_miRNA.close()

    # make a list to the store the UTR length of miRNA targets 
    UTR_length = []
    
    # create lists to store the number of miRNAs, binding sites and sites per miRNA for the targets
    N_mirnas = []
    N_sites = []
    ratio = []
    
    # populate the different lists keeping the same order within each list
    for gene in targets:
            UTR_length.append(CDS_UTR[gene])
            N_mirnas.append(targets[gene][0])
            N_sites.append(targets[gene][1])
            ratio.append(targets[gene][2])
                    
    return [UTR_length, N_mirnas, N_sites, ratio]

    



def correlation_CDSUTR_length_miRNAs(targetscan_CDS_UTR_input_file, species, CNV_miRNA_file):
    '''
    (file, str, file) -> tuple
    Performs a Spearman's rank correlation between the UTR length of each target and:
    number the miRNAs, number of sites and ratio of sites per miRNA of the same target
    Returns a tuple of tuple with each inner tuple containing the rho and p values
    '''

    # get the list of UTR length, and N mirnas, sites and ratio
    UTR_length, N_mirnas, N_sites, ratio = make_CDSUTR_length_miRNA_regulation_lists(targetscan_CDS_UTR_input_file, species, CNV_miRNA_file)

    # perform Sperman's rank correlation between UTR length and number of miRNAs
    rho_mirnas, p_mirnas = stats.spearmanr(UTR_length, N_mirnas)

    # perform Sperman's rank correlation between UTR length and number of sites
    rho_sites, p_sites = stats.spearmanr(UTR_length, N_sites)

    # perform Sperman's rank correlation between UTR length and ratio of sites per mirna
    rho_ratio, p_ratio = stats.spearmanr(UTR_length, ratio)

    return (rho_mirnas, p_mirnas), (rho_sites, p_sites), (rho_ratio, p_ratio)



def compare_CDSUTR_length_CNV_nonCNV(targetscan_CDS_UTR_input_file, species, CNV_miRNA_file):
    '''
    (file, str, file) -> tuple of tuples
    Returns a tuple of tuples with the mean CDS+UTR length and SEM for CNV genes, and a tuple for non-CNV genes
    and a tuple with the Wilcoxon ranksum z-score and its p-value for testing
    the length difference of CDS + UTR between CNV and nonCNV genes
    '''

    # get the sequence of all CDS + UTRs for the given species
    CDS_UTR = {}

    # open file for reading
    cds = open(targetscan_CDS_UTR_input_file, 'r')
    for line in cds:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            seq = line[-1]
            CDS_UTR[gene] = len(seq)
    #close file after reading
    cds.close()

    # make a dictionnary with information for each target transcript
    targets = {}
    CNV_miRNA = open(CNV_miRNA_file, 'r')
    CNV_miRNA.readline()
    for line in CNV_miRNA:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if species == 'fly':
                # use the transcript as key
                targets[line[1]] = line[-1]
            else:
                # use the gene as key
                targets[line[0]] = line[-1]
    # close after reading
    CNV_miRNA.close()

    # make a list to the store the CDS + UTR length of CNV and nonCNV genes
    CNV_length = []
    nonCNV_length = []
    
    # populate the different lists keeping the same order within each list
    for gene in targets:
        if targets[gene] == 'CNV':
            CNV_length.append(CDS_UTR[gene])
        elif targets[gene] == 'non-CNV':
            nonCNV_length.append(CDS_UTR[gene])

    # perform a Wilcoxon rank sum test
    diff = stats.ranksums(CNV_length, nonCNV_length)

    # compute means and SEM
    mean_CNV = stats.tmean(CNV_length)
    sem_CNV = stats.sem(CNV_length)
    mean_nonCNV = stats.tmean(nonCNV_length)
    sem_nonCNV = stats.sem(nonCNV_length)
                    
    return (mean_CNV, sem_CNV), (mean_nonCNV, sem_nonCNV), diff
