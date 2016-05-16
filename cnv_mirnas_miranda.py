def transcripts_gene_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score):
    '''
    (file, file, file, file, str) -> dict
    Returns a dictionnary with the transcript name as key and gene name as value
    with option to select only good_score sites or all sites including poor_score sites
    '''

    # make a dictionnary to store the transcript : gene pairs
    transcripts_genes = {}

    # open and read the files    
    files = [conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor]
    if site_score == 'all_sites':
        open_files = files[0:]
    elif site_score == 'good_scores':
        open_files = files[0:2]

    for file in open_files:
        targets = open(file, 'r')
        header = targets.readline()
        for line in targets:
            if line.rstrip() != '':
                line = line.split()
                transcript = line[4]
                gene = line[3]
                transcripts_genes[transcript] = gene
        targets.close()
    
    return transcripts_genes
    

def miRNA_regulation_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score):
    '''
    (file, file, file, file, str) -> dict
    Returns a dictionnary of dictionnaries for each transcript with the the inner dictionnaries containing the number of sites for each mirna
    with option to select only good_score sites or all sites including poor_score sites
    '''

    # make a dictionnary to store the transcripts {transcript1: {mirna1: N_sites, mirna2: N_sites}}
    regulated_transcripts = {}

    # open and read the files    
    files = [conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor]
    if site_score == 'all_sites':
        open_files = files[0:]
    elif site_score == 'good_scores':
        open_files = files[0:2]

    for file in open_files:
        targets = open(file, 'r')
        header = targets.readline()
        for line in targets:
            if line.rstrip() != '':
                line = line.split()
                transcript = line[4]
                mirna = line[1]
                if transcript in regulated_transcripts:
                    if mirna in regulated_transcripts[transcript]:
                        regulated_transcripts[transcript][mirna] += 1
                    else:
                        regulated_transcripts[transcript][mirna] = 1
                else:
                    regulated_transcripts[transcript] = {}
                    regulated_transcripts[transcript][mirna] = 1
        targets.close()
       
    return regulated_transcripts


def worm_make_miRNA_regulators_table_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score,
                                             non_conconserved_poor, site_score, all_gene_file, CNV_file, outputfile):
    '''
    (file, file, file, file, str, file, file, file) -> file
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites,
    the mean number of sites / mirna and whether the gene is in a CNV or not, with option to include only good_scores sites or all_sites
    '''
    # get the transcript : gene pairs
    transcripts_genes = transcripts_gene_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)

    # get valid genes
    valid_genes = get_valid_worm_gene_ID(all_gene_file)

    # get the CNV genes
    CNV_genes = get_worm_CNV_genes(CNV_file, all_gene_file)

    # get the number of regulators for each transcript
    regulated_transcripts = miRNA_regulation_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)

    # make a dictionnary {gene1:[N_mirnas, N_sites, ratio, transcript, CNV]
    regulated_genes = {}
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
        for wormbase_gene_ID in valid_genes:
            if gene in valid_genes[wormbase_gene_ID]:
                if wormbase_gene_ID in CNV_genes:
                    regulated_genes[gene].append('CNV')
                else:
                    regulated_genes[gene].append('non-CNV')
                break

    # make a set of valid gene names and valid gene symbols
    valid_gene_names = set()
    for wormbase_gene_ID in valid_genes:
        for gene in valid_genes[wormbase_gene_ID]:
            valid_gene_names.add(gene)
    
    # remove target genes that are not valid
    to_remove = []
    for gene in regulated_genes:
        if gene not in valid_gene_names:
            to_remove.append(gene)
    if len(to_remove) != '':
        for gene in to_remove:
            del regulated_genes[gene]
    
    # write to file
    newfile = open(outputfile, 'w')
    newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
    for gene in regulated_genes:
        newfile.write(gene + '\t' + regulated_genes[gene][-2] + '\t' + str(regulated_genes[gene][0]) + '\t' + str(regulated_genes[gene][1]) + '\t'
                      + str(regulated_genes[gene][2]) + '\t' + regulated_genes[gene][-1] + '\n')

    newfile.close()

    

def get_valid_worm_gene_ID(all_gene_file):
    '''
    (file) -> dict
    Return a dictionnary of valid C. elegans genes that are alive in the current genome annotation file
    with WBGene ID as key and a list of gene ID and gene name when name is available
    '''

    # make a dictionnary {WBGene: [gene_ID, gene_name]}
    valid_genes = {}
    IDs = open(all_gene_file, 'r')
    for line in IDs:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            if line[-1] == 'Live':
                while '' in line:
                    line.remove('')
                if len(line) == 4:
                    valid_genes[line[1]] = [line[2]]
                elif len(line) == 5:
                    valid_genes[line[1]] = [line[2], line[3]]

    IDs.close()
    return valid_genes


def get_worm_CNV_genes(CNV_file, all_gene_file):
    '''
    (file, file) -> set
    Returns a set of CNV genes using the WBGene ID
    '''
    
    # get the dictionnary of valid worm genes
    valid_genes = get_valid_worm_gene_ID(all_gene_file)

    # make a set of CNV genes using the WBGene ID
    CNV = open(CNV_file, 'r')
    CNV_genes = set()
    header = CNV.readline()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            line = line.split()
            CNV_genes.add(line[2])

    # remove non valid CNV genes
    to_remove = []
    for gene in CNV_genes:
        if gene not in valid_genes:
            to_remove.append(gene)
    for gene in to_remove:
        CNV_genes.discard(gene)

    CNV.close()
    return CNV_genes


def human_make_miRNA_regulators_table_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score,
                                             non_conconserved_poor, site_score, all_gene_file, CNV_file, outputfile):
    '''
    (file, file, file, file, str, file, file, file) -> file
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites,
    the mean number of sites / mirna and whether the gene is in a CNV or not, with option to include only good_scores sites or all_sites
    '''

    from CNV_miRNAs_TargetScan import sort_valid_human_genes
    from CNV_miRNAs_TargetScan import human_CNV_genes
    
    # get the transcript : gene pairs
    transcripts_genes = transcripts_gene_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)

    # get valid genes
    valid_genes = sort_valid_human_genes(all_gene_file)

    # get the CNV genes
    CNV_genes = human_CNV_genes(CNV_file, all_gene_file)

    # get the number of regulators for each transcript
    regulated_transcripts = miRNA_regulation_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)

    # make a dictionnary {gene1:[N_mirnas, N_sites, ratio, transcript, CNV]
    regulated_genes = {}
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
    
    # write to file
    newfile = open(outputfile, 'w')
    newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
    for gene in regulated_genes:
        newfile.write(gene + '\t' + regulated_genes[gene][-2] + '\t' + str(regulated_genes[gene][0]) + '\t' + str(regulated_genes[gene][1]) + '\t'
                      + str(regulated_genes[gene][2]) + '\t' + regulated_genes[gene][-1] + '\n')

    newfile.close()


def fly_make_miRNA_regulators_table_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score,
                                             non_conconserved_poor, site_score, fly_gene_coordinates, CNV_file, outputfile):
    '''
    (file, file, file, file, str, file, file, file) -> file
    For a single transcript / gene, write the number of miRNA regulators, number of miRNA binding sites,
    the mean number of sites / mirna and whether the gene is in a CNV or not, with option to include only good_scores sites or all_sites
    '''

    from CNV_miRNAs_TargetScan import sort_valid_fly_genes
        
    # get the transcript : gene pairs
    transcripts_genes = transcripts_gene_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)

    # get valid genes
    valid_genes = sort_valid_fly_genes(fly_gene_coordinates)

    # get the CNV genes
    CNV_genes = set()
    CNV = open(CNV_file, 'r')
    CNV.readline()
    for line in CNV:
        line = line.rstrip()
        if line != '':
            line = line.split()
            CNV_genes.add(line[0])
            CNV_genes.add(line[1])
    
    # get the number of regulators for each transcript
    regulated_transcripts = miRNA_regulation_miranda(conserved_good_score, non_conserved_good_score, conserved_poor_score, non_conconserved_poor, site_score)

    # make a dictionnary {gene1:[N_mirnas, N_sites, ratio, transcript, CNV]
    regulated_genes = {}
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
    
    # write to file
    newfile = open(outputfile, 'w')
    newfile.write('gene' + '\t' + 'transcript' + '\t' + 'N_mirnas' + '\t' + 'N_sites' + '\t' + 'N_sites_per_miRNA' + '\t' + 'CNV' + '\n')
    for gene in regulated_genes:
        newfile.write(gene + '\t' + regulated_genes[gene][-2] + '\t' + str(regulated_genes[gene][0]) + '\t' + str(regulated_genes[gene][1]) + '\t'
                      + str(regulated_genes[gene][2]) + '\t' + regulated_genes[gene][-1] + '\n')
    CNV.close()
    newfile.close()
