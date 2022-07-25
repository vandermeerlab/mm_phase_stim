# Script to prefix .t, .ntt and ExpKeys files with the correct Suffix

import os
import shutil

# Change this to the folder where you want the renaming to happen
target_dir = "E:/Dropbox (Dartmouth College)/manish_data/M321/M321-2022-07-13"
# Might need to hard-code this in some-cases 
session_name = target_dir.split('/')[-1]
os.chdir(target_dir)
all_files = os.listdir()
all_files = [file for file in all_files if os.path.isfile(file) and \
 ('.t' in file or '.ntt' in file)]
for file in all_files:
    new_filename = '-'.join([session_name, file])
    shutil.move(file, new_filename)   
