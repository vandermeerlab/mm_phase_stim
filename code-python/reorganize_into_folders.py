# Script to organize session files into folders

import os
import shutil

#Change this to the folder where you want the reorganization to happen
target_dir = "D:/Dropbox (Dartmouth College)/NACC_BCOD/Cohort_2/M353"
os.chdir(target_dir)
all_files = os.listdir()
all_files = [file for file in all_files if os.path.isfile(file)]
for file in all_files:
    toks = file.split('-')
    folder = '-'.join(toks[0:-1])
    if not os.path.isdir(folder):
        os.mkdir(folder)
    shutil.move(file, '/'.join([folder,file]))   