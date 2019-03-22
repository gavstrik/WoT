import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'

"""
This python script plotting summary statistics for
the file-condition.
"""

datafiles = [
            '../data/dots.xls',
            ]


# outliers with an error rate above 10 are removed
def remove_outliers(data, true_number_of_dots):
    max_error_rate = 10
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    newdata = [i for i in data if i >= lower and i <= upper]
    return np.array(newdata)


def plot_stats(df_all):
    colors = ['#fdae6b', '#f16913', '#d94801', '#7f2704']

    # global figure settings:
    fig, axarr = plt.subplots(nrows=1, figsize=(11.4,2.5))
    plt.xlabel('log(estimate)', fontsize=18)

    # plot
    lv = []
    lv.append('')
    df_control = df_all[df_all.v == 0]
    control_guesses = remove_outliers(df_control.guess.values, 1097)
    true_number_of_dots = 1097
    for view_idx, m in enumerate([(0, 'history'), (9, 'file'), (9, 'history')]):
        df = df_all[(df_all.v == m[0]) & (df_all.method == m[1])]
        guesses = remove_outliers(df.guess.values, true_number_of_dots)
        if m == (9, 'file'):
            g = remove_outliers(df.guess.values, true_number_of_dots)

        # calculate the p-value using the two-sided Mann-Whitney U test
        _, pvalue = mannwhitneyu(control_guesses, guesses)
        if 0.01 < pvalue <= 0.05:
            lv.append(str(m[1])+', v='+str(m[0])+' (N='+str(len(guesses))+')*')
        elif 0.001 < pvalue <= 0.01:
            lv.append(str(m[1])+', v='+str(m[0])+' (N='+str(len(guesses))+')**')
        elif pvalue <= 0.001:
            lv.append(str(m[1])+', v='+str(m[0])+' (N='+str(len(guesses))+')***')
        else:
            lv.append(str(m[1])+', v='+str(m[0])+' (N='+str(len(guesses))+')')
        print(m[0], m[1], pvalue)
        if m == (9, 'history'):
            print('comparing file with history:', mannwhitneyu(guesses, g))

        # calculate quartiles and deciles
        quan_min = np.log(np.quantile(guesses, 0.25))
        quan_max = np.log(np.quantile(guesses, 0.75))
        dec_min = np.log(np.quantile(guesses, 0.1))
        dec_max = np.log(np.quantile(guesses, 0.9))

        # plot the thread stats for the given number of views
        axarr.plot(np.log(np.median(guesses)),
                             [view_idx], 'o', markersize=9, c=colors[view_idx])
        axarr.plot(np.log(np.mean(guesses)),
                             [view_idx], 'v', markersize=7, c=colors[view_idx])
        axarr.hlines(y=view_idx, xmin=quan_min, xmax=quan_max,
                               linewidth=4, color=colors[view_idx])
        axarr.hlines(y=view_idx, xmin=dec_min, xmax=dec_max,
                               linewidth=2, color=colors[view_idx])

    # plotting paraphernalia for the subfigure for the number of dots
    axarr.set_yticks([-0.5, 0, 1, 2, 2.5])
    axarr.tick_params(axis="y", length=0)  # remove small ticks
    axarr.set_yticklabels(lv, fontsize=14, ha='left')
    axarr.tick_params(axis='y', direction='out', pad=175)
    axarr.axvline(x=np.log(true_number_of_dots),
                               linewidth=2, color='#34495e')
    axr = axarr.twinx()
    axr.set_ylabel('d='+str(true_number_of_dots), fontsize=18)
    plt.yticks([])  # removing axis numbers with empty list
    axarr.tick_params(axis="x", length=0)
    axarr.grid(True)
    axarr.set_xticklabels(['','5.5','6.0','6.5','7.0','7.5','8',''], fontdict={'fontsize': 14})
    axarr.set_xlim([5.5, 8])  # set x-axis dimensions
    plt.tight_layout()

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    plt.savefig(os.path.join(PLOTS_DIR, 'figS10.png'), transparent=True, dpi=300)
    plt.savefig(os.path.join(PLOTS_DIR, 'figS10.pdf'), transparent=True, dpi=300)
    plt.show()


# load the data
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_excel(datafile))
    df_all = df_all.append(df, sort=True)
df_all = df_all[(df_all.session == '699h8rze') | (df_all.session == 'hal5jdl0') | (df_all.session == 'wv4xujg7')]

plot_stats(df_all)
