# Script to reorganize split_ntt_files into appropriate folders

import os
import shutil

#Change this to the folder where you want the reorganization to happen
target_dir = "D:/Dropbox (Dartmouth College)/NACC_BCOD/Cohort_2/M333/M333-2022-05-27"
os.chdir(target_dir)
all_files = os.listdir()
all_files = [file for file in all_files if os.path.isfile(file) and 'split' in file]
for file in all_files:
    toks = file.split('_')
    folder = '-'.join(['Recording', toks[-1][0]])
    if not os.path.isdir(folder):
        os.mkdir(folder)
    new_filename = toks[-2] + '.ntt'  
    shutil.move(file, '/'.join([folder,new_filename]))   