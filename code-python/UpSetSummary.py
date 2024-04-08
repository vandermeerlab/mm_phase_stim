import upsetplot as up
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os

# This is to ensure that the text in the figure exported by matplotlib is editable in Adobe Illustrator
mpl.rcParams['pdf.fonttype'] = 42 

# Change the font to Arial 
mpl.rcParams['font.family'] = 'Arial'

# Change the font size to 30
mpl.rcParams['font.size'] = 30


filepath1 = os.path.join('C:','Users','mvdmlab','Desktop','dStr_sig.csv')
filepath2 = os.path.join('C:','Users','mvdmlab','Desktop','vStr_sig.csv')
filepath3 = os.path.join('C:','Users','mvdmlab','Desktop','all_sig.csv')

dStr_data = pd.read_csv(filepath1, index_col=[0,1,2,3])
vStr_data = pd.read_csv(filepath2, index_col=[0,1,2,3])
all_data = pd.read_csv(filepath3, index_col=[0,1,2,3])


fig = plt.figure(figsize=(21, 15))
ax = up.plot(dStr_data, fig=fig, sum_over='Count', sort_by='input', sort_categories_by='input', show_counts=True, show_percentages=True, orientation = 'vertical', element_size=None)
ax['totals'].set_yticks([0,15])
ax['totals'].set_ylabel('Totals')
ax['intersections'].set_ylabel('Count')
ax['intersections'].set_xticks([0,15])
ax['totals'].tick_params('both', length=20, which='major')
ax['intersections'].tick_params('both', length=20, which='major')
plt.suptitle('dStr')
# Save figure
plt.savefig(os.path.join('C:','Users','mvdmlab','Desktop','dStr_summary.pdf'))

fig = plt.figure(figsize=(21, 15))
ax = up.plot(vStr_data, fig=fig, sum_over='Count', sort_by='input', sort_categories_by='input', show_counts=True, show_percentages=True, orientation = 'vertical', element_size=None)
ax['totals'].set_yticks([0,36])
ax['totals'].set_ylabel('Totals')
ax['intersections'].set_ylabel('Count')
ax['intersections'].set_xticks([0,36])
ax['totals'].tick_params('both', length=20, which='major')
ax['intersections'].tick_params('both', length=20, which='major')
plt.suptitle('vStr')
# Save figure
plt.savefig(os.path.join('C:','Users','mvdmlab','Desktop','vStr_summary.pdf'))

fig = plt.figure(figsize=(21, 15))
ax = up.plot(all_data, fig=fig, sum_over='Count', sort_by='input', sort_categories_by='input', show_counts=True, show_percentages=True, orientation = 'vertical', element_size=None)
ax['totals'].set_yticks([0,51])
ax['totals'].set_ylabel('Totals')
ax['intersections'].set_ylabel('Count')
ax['intersections'].set_xticks([0,51])
ax['totals'].tick_params('both', length=20, which='major')
ax['intersections'].tick_params('both', length=20, which='major')
# Save figure
plt.savefig(os.path.join('C:','Users','mvdmlab','Desktop','all_summary.pdf'))