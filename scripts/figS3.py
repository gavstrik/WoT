import numpy as np
import pandas as pd
from statsmodels.graphics.gofplots import qqplot
import matplotlib.pyplot as plt
import matplotlib
# plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"

"""
QQ-plots.
"""

datafiles = [
            '../data/dots/all_dots_trimmed_anonymous.csv',
            '../data/ox/all_ox_trimmed_anonymous.csv',
            ]


def basic_stats(df1):
    fig, axes = plt.subplots(nrows=3, ncols=2,sharex=True, figsize=(11, 11))
    axes = axes.flatten()
    images = [(55, 3), (148, 9), (403, 3), (1097, 0), (1233, 0), (1233, 9)]
    for image_idx, (d, views) in enumerate(images):
        df = df1[(df1['dots'] == d) & (df1['views'] == views)]
        guesses = df['guess'].values
        X = np.log(guesses)
        qqplot(X, line='s', ax=axes[image_idx])
        axes[image_idx].set_title('d = '+str(d)+', v = '+str(views)+', N ='+str(len(guesses)))
        axes[image_idx].set_xlabel('')
    plt.tight_layout()
    # Remember: save as pdf and transparent=True for Adobe Illustrator
    fig.savefig('../plots/FigS3.pdf', transparent=True, dpi=300)
    plt.show()

# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)
df_all = df_all[df_all['method'] == 'history']

basic_stats(df_all)
