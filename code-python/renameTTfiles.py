# Script to convert TTX* files into TT0X* files

import os
import shutil

#Change this to the folder where you want the renaming to happen
target_dir = "E:/Dropbox (Dartmouth College)/manish_data/M321/M321-2022-07-08"
os.chdir(target_dir)
all_files = os.listdir()
all_files = [file for file in all_files if os.path.isfile(file) and 'TT' in file and 'TT0' not in file]
for file in all_files:
    # Skip file if TT number greater than 9 (TT10 etc)
    toks = file.split('TT')
    if len(toks[1]) > 1 and toks[1][1].isnumeric():
        continue
    new_filename = ''.join(['TT0', toks[-1]])
    shutil.move(file, new_filename)   
