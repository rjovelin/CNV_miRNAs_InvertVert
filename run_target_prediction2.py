import os



# run prediction for worm
os.system('perl targetscan_60.pl Celegans_mirna_input.txt Celegans_CDS_UTR_input.txt Celegans_target_sites_CDS_UTR.txt')

# run prediction for drosophila
os.system('perl targetscan_60.pl Drosophila_mirna_input.txt Drosophila_CDS_UTR_input.txt Drosophila_target_sites_CDS_UTR.txt')

# run prediction for zebrafish
os.system('perl targetscan_60.pl Zebrafish_mirna_input.txt Zebrafish_CDS_UTR_input.txt Zebrafish_target_sites_CDS_UTR.txt')

# run prediction for human
os.system('perl targetscan_60.pl human_mirna_input.txt Human_CDS_UTR_input.txt human_target_sites_CDS_UTR.txt')





