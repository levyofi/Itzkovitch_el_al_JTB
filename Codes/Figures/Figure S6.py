import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

mod1 = 1 # model 1 == before ML
mod2 = 3 # model 3 == after ML + mean(PE)

df_map = pd.read_csv(f'../Example data/Tables/SampledPixels_13.csv')

fig, ax = plt.subplots(2,3, figsize = (18,10))
colors = ['lightblue', '#FFD580', '#CF9FFF']
c1 = colors[mod1-1]
c2 = colors[mod2-1]

names = ['before ML', 'after ML', 'after ML + mPE']
name1 = names[mod1-1]
name2 = names[mod2-1]
# Physical prediction
ax11 = ax[0,0] #plt.subplot2grid((7, 14), (0, 0), colspan=4, rowspan=3)
sns.scatterplot(data = df_map, x = 'Phy', y = 'M1', s = 1, color = c1, label = 'Physical', ax = ax11, legend = False)
sns.scatterplot(data = df_map, x = 'Phy', y = 'M2', s = 1, color = c2, label = 'Statistical', ax = ax11, legend = False)
ax11.set_ylabel(u'Prediction Error (\u00B0C)', size = 16)
ax11.set_xlabel(u'Predicted Ground Temperature (K)', size = 16)
ax11.hlines(0, 292, 330, color = 'grey', linestyle = 'dashed')

# TGI
ax21 = ax[0,1] #plt.subplot2grid((7, 14), (0, 5), colspan=4, rowspan=3)
sns.scatterplot(data = df_map, x = 'TGI', y = 'M1', s = 1, color = c1, label = 'Physical', ax = ax21, legend = False)
sns.scatterplot(data = df_map, x = 'TGI', y = 'M2', s = 1, color = c2, label = 'Statistical', ax = ax21, legend = False)
ax21.set_ylabel('')
ax21.set_xlabel('TGI', size = 16)
ax21.set_xticks([-0.1, -0.05, 0, 0.04])
ax21.hlines(0, -0.1, 0.04, color = 'grey', linestyle = 'dashed')
ax21.set_xlim(-0.1,0.04)

# Height
ax31 = ax[0,2] #plt.subplot2grid((7, 14), (0, 10), colspan=4, rowspan=3)
sns.scatterplot(data = df_map, x = 'Height', y = 'M1', s = 1, color = c1, label = 'Physical', ax = ax31, legend = False)
sns.scatterplot(data = df_map, x = 'Height', y = 'M2', s = 1, color = c2, label = 'Statistical', ax = ax31, legend = False)
ax31.set_ylabel('')
ax31.set_xlabel('Height (m)', size = 16)
ax31.hlines(0, 0, 30, color = 'grey', linestyle = 'dashed')

# Skyview
ax41 = ax[1,0] #plt.subplot2grid((7, 14), (4, 0), colspan=4, rowspan=3)
sns.scatterplot(data = df_map, x = 'Skyview', y = 'M1', s = 1, color = c1, label = 'Physical', ax = ax41, legend = False)
sns.scatterplot(data = df_map, x = 'Skyview', y = 'M2', s = 1, color = c2, label = 'Statistical', ax = ax41, legend = False)
ax41.set_ylabel(u'Prediction Error (\u00B0C)', size = 16)
ax41.set_xlabel('Skyview (dec %)', size = 16)
ax41.hlines(0, 0, 1, color = 'grey', linestyle = 'dashed')

# Shade
ax51 = ax[1, 1]  # Reference subplot
sns.boxplot(
    data=df_map.melt('Shade', value_vars=['M1', 'M2'], var_name='Model'),
    x='Shade', y='value', hue='Model',
    ax=ax51, palette={"M1": c1, "M2": c2}, showfliers=False  # Disable outliers (optional)
)
# Adjust labels and legend
ax51.set_xlabel('Shade (yes/no)', size=16)
ax51.hlines(0, -0.5, 1.5, color='grey', linestyle='dashed')  # Adjust x-limits for hline
ax51.set_xticklabels(['no', 'yes'])
ax51.set_ylabel('')

# Remove duplicate legend
handles, labels = ax51.get_legend_handles_labels()
ax51.legend(handles[:2], [name1, name2], loc='lower center', ncol=2, bbox_to_anchor=(0.5, -0.35))

# Real Solar
ax61 = ax[1,2] #plt.subplot2grid((7, 14), (4, 10), colspan=4, rowspan=3)
sns.scatterplot(data = df_map, x = 'RealSolar', y = 'M1', s = 1, color = c1, label = 'Physical', ax = ax61, legend = False)
sns.scatterplot(data = df_map, x = 'RealSolar', y = 'M2', s = 1, color = c2, label = 'Statistical', ax = ax61, legend = False)
ax61.set_xlabel(r'Solar Radiation (Wm$^{-2}$)', size=16)
ax61.set_ylabel('')
ax61.hlines(0, 0, 1000, color = 'grey', linestyle = 'dashed')

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.2)
fig.savefig('Figure S6.png', dpi=300, bbox_inches='tight')
