import upsetplot as up
import matplotlib.pyplot as plt
import pandas as pd
import os

filepath1 = os.path.join('C:','Users','mvdmlab','Desktop','dStr_sig.csv')
filepath2 = os.path.join('C:','Users','mvdmlab','Desktop','vStr_sig.csv')

dStr_data = pd.read_csv(filepath1, index_col=[0,1,2])
vStr_data = pd.read_csv(filepath2, index_col=[0,1,2])
fig = plt.figure(figsize=(30, 10))
ax = up.plot(dStr_data, fig=fig, sum_over='Count', sort_by='input', sort_categories_by='-input')
ax['totals'].set_visible(False)
ax['intersections'].set_ylabel('Count')
plt.suptitle('dStr')
plt.savefig(os.path.join('C:','Users','mvdmlab','Desktop','dStr_summary.svg'))

fig = plt.figure(figsize=(30, 10))
ax = up.plot(vStr_data, fig=fig, sum_over='Count', sort_by='input', sort_categories_by='-input')
ax['totals'].set_visible(False)
ax['intersections'].set_ylabel('Count')
plt.suptitle('vStr')
plt.savefig(os.path.join('C:','Users','mvdmlab','Desktop','vStr_summary.svg'))