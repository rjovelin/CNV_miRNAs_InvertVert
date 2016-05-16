from CNV_miRNAs_TargetScan import *


def make_targetscan_UTR_input(CNV_miRNA_file, UTR_file, species, threshold, outputfile):
    '''
    (file, file, str, str, file) -> file
    Get the UTR sequence for all the genes in the CNV_miRNA file, remove gaps and
    chop all sequence to 1 Kb if length > 1 Kb and do not keep if length < 1 Kb
    Save sequences with same length to file as input for running targetscan
    '''

    # open miRNAs CNV file for reading
    CNVs = open(CNV_miRNA_file, 'r')
    # skip header
    CNVs.readline()

    # make a dictionnary of {transcript : [gene, CNV_status]}
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
    # close file after reading
    CNVs. close()

    # get the sequence of all UTRs for the given species {transcript : UTR}
    all_UTR = UTR_sequence(UTR_file, species)

    # create a new dictionnary to hold the chopped UTR of same length
    chopped_UTR = {}
    # populate the new dict
    for transcript in CNV_genes:
        sequence = all_UTR[transcript]
        sequence = sequence.replace('-', '').upper()
        if len(sequence) > threshold:
            sequence = sequence[:threshold]
            chopped_UTR[transcript] = sequence

    # open outputfile for writing
    newfile = open(outputfile, 'w')
    for transcript in chopped_UTR:
        newfile.write(CNV_genes[transcript][0] + '\t')
        if species == 'worm':
            newfile.write('6239' + '\t')
        elif species == 'fly':
            newfile.write('7227' + '\t')
        elif species == 'human':
            newfile.write('9606' + '\t')
        elif species == 'zebrafish':
            newfile.write('7955' + '\t')
        newfile.write(chopped_UTR[transcript] + '\n')
    # close file after writing
    newfile.close()
        




def mirna_family_seed_name(miRNAfamily_info, species):
    '''
    (file, str) -> dict
    Return a dictionnary with the seed as key and a list as value containing the name of the family and the species ID
    '''

    # open the mirfam info file for reading
    mirna = open(miRNAfamily_info, 'r')
    mirna.readline()

    # create a dictionnary to hold the seed and family name
    family = {}

    # go through the file and populate the dictionnary
    for line in mirna:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            seed = line[1]
            mirfam = line[0]
            if species == 'human':
                if line[2] == '9606':
                    family[seed] = [mirfam, '9606']
            elif species == 'worm':
                if line[2] == '6239':
                    family[seed] = [mirfam, '6239']
            elif species == 'fly':
                family[seed] = [mirfam, '7227']
            elif species == 'zebrafish':
                if line[2] == '7955':
                    family[seed] = [mirfam, '7955']
    # close file after reading
    mirna.close()

    return family




def mirnas_per_family(miRNAfamily_info, species):
    '''
    (file, str) -> dict
    Return a dictionnary with the seed as key and set of mirnas in the same family
    '''

    # open the mirfam info file for reading
    mirna_file = open(miRNAfamily_info, 'r')
    mirna_file.readline()

    # create a dictionnary to hold the seed and family name
    same_family = {}

    # go through the file and populate the dictionnary
    for line in mirna_file:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # grab the seed and the mirna
            seed = line[1]
            mirna = line[3]
            # populate the dictionnary if correct species
            if species == 'human':
                if line[2] == '9606':
                    # add the mirna to the set if seed already in dict
                    if seed in same_family:
                        same_family[seed].add(mirna)
                    # if seed not in dict, add the seed : set(mirna) pair
                    else:
                        same_family[seed] = {mirna,}
            elif species == 'worm':
                if line[2] == '6239':
                    # add the mirna to the set if seed already in dict
                    if seed in same_family:
                        same_family[seed].add(mirna)
                    # if seed not in dict, add the seed : set(mirna) pair
                    else:
                        same_family[seed] = {mirna,}
            elif species == 'fly':
                mirna = line[2]
                # add the mirna to the set if seed already in dict
                if seed in same_family:
                    same_family[seed].add(mirna)
                # if seed not in dict, add the seed : set(mirna) pair
                else:
                    same_family[seed] = {mirna,}
            elif species == 'zebrafish':
                if line[2] == '7955':
                    # add the mirna to the set if seed already in dict
                    if seed in same_family:
                        same_family[seed].add(mirna)
                    # if seed not in dict, add the seed : set(mirna) pair
                    else:
                        same_family[seed] = {mirna,}
 
    # close file after reading
    mirna_file.close()

    return same_family

    

# make mirna file as input for targetscan    

def make_targetscan_mirna_imput(miRNAfamily_info, species, outputfile):
    '''
    (file, str, file) -> file
    Save the miRNA family info as input for running Targetscan
    '''

    # get the mirna family names and seeds
    family = mirna_family_seed_name(miRNAfamily_info, species)

    # open outputfile for writing
    newfile = open(outputfile, 'w')
    for seed in family:
        newfile.write(family[seed][0] + '\t' + seed + '\t' + family[seed][1] + '\n')
    # close file after writing
    newfile.close()
        

def make_summary_table_CNV_miRNA(targetscan_output, CNV_miRNA_file, miRNAfamily_info, species, outputfile):
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

    # make a dictionnary with the gene name as key and a list with the transcript name and CNV status as value
    CNV_genes = {}
    
    # open miRNAs CNV file for reading
    CNVs = open(CNV_miRNA_file, 'r')
    # skip header
    CNVs.readline()

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
        newfile.write(gene + '\t' + CNV_genes[gene][0] + '\t' + str(mirna_count) + '\t' + str(site_count) + '\t')
        if mirna_count != 0:
            newfile.write(str(site_count / mirna_count) + '\t')
        else:
            newfile.write('NA' + '\t')
        newfile.write(CNV_genes[gene][1] + '\n')
    # close file ater writing
    newfile.close()
        
            
            

    
            

    
