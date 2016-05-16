import os


# run prediction for worm
os.system('perl targetscan_60.pl Celegans_mirna_input.txt Celegans_300bp_UTR_input.txt Celegans_target_sites_300bp_UTR.txt')
os.system('perl targetscan_60.pl Celegans_mirna_input.txt Celegans_500bp_UTR_input.txt Celegans_target_sites_500bp_UTR.txt')
os.system('perl targetscan_60.pl Celegans_mirna_input.txt Celegans_700bp_UTR_input.txt Celegans_target_sites_700bp_UTR.txt')

# run prediction for drosophila
os.system('perl targetscan_60.pl Drosophila_mirna_input.txt Drosophila_300bp_UTR_input.txt Drosophila_target_sites_300bp_UTR.txt')
os.system('perl targetscan_60.pl Drosophila_mirna_input.txt Drosophila_500bp_UTR_input.txt Drosophila_target_sites_500bp_UTR.txt')
os.system('perl targetscan_60.pl Drosophila_mirna_input.txt Drosophila_700bp_UTR_input.txt Drosophila_target_sites_700bp_UTR.txt')

# run prediction for zebrafish
os.system('perl targetscan_60.pl Zebrafish_mirna_input.txt Zebrafish_300bp_UTR_input.txt Zebrafish_target_sites_300bp_UTR.txt')
os.system('perl targetscan_60.pl Zebrafish_mirna_input.txt Zebrafish_500bp_UTR_input.txt Zebrafish_target_sites_500bp_UTR.txt')
os.system('perl targetscan_60.pl Zebrafish_mirna_input.txt Zebrafish_700bp_UTR_input.txt Zebrafish_target_sites_700bp_UTR.txt')

# run prediction for human
os.system('perl targetscan_60.pl human_mirna_input.txt human_300bp_UTR_input.txt human_target_sites_300bp_UTR.txt')
os.system('perl targetscan_60.pl human_mirna_input.txt human_500bp_UTR_input.txt human_target_sites_500bp_UTR.txt')
os.system('perl targetscan_60.pl human_mirna_input.txt human_700bp_UTR_input.txt human_target_sites_700bp_UTR.txt')
