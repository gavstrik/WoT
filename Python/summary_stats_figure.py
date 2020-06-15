import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'

"""
This python script plotting summary statistics for
the WoT-experiments on Amazon Mechanical Turk.
"""

datafiles = [
            '../CleanedData/dots.xlsx',
            ]


def plot_stats(df_all):
    dots = [1097, 403, 148, 55]
    views = [0, 1, 3, 9]
    # colors = ['#fdae6b', '#f16913', '#d94801', '#7f2704']
    color_greys = ['#252525', '#525252', '#737373', '#969696', '#bdbdbd']
    color_blues = ['#08519c', '#2171b5', '#4292c6', '#6baed6', '#9ecae1']
    color_reds = ['#a50f15', '#cb181d', '#ef3b2c', '#fb6a4a','#fc9272']

    # global figure settings:
    fig, axarr = plt.subplots(nrows=4, ncols=2, sharex=True, figsize=(15.6, 8.45))
    # plt.rcParams.update({"yticklabel.direction" : "out"})

    # plot
    for idm, method in enumerate(['history', 'max']):
        dfm = df_all[(df_all.method == method) | (df_all.v == 0)]
        lv = [[''] for i in range(5)]
        for true_dots_idx, true_dots in enumerate(dots):
            df_i = dfm[dfm.d == true_dots]
            df_control = df_i[df_i.v == 0]
            # control_guesses = remove_outliers(df_control.guess.values, true_dots)
            control_guesses = df_control.guess.values
            true_number_of_dots = true_dots
            for view_idx, view in enumerate(views):
                df = df_i[df_i.v == view]
                # guesses = remove_outliers(df.guess.values, true_dots)
                guesses = df.guess.values
                if view_idx == 0:
                    colors = color_greys
                else:
                    if idm == 0:
                        colors = color_blues
                    else:
                        colors = color_reds

                # calculate the p-value using the two-sided Mann-Whitney U test
                _, pvalue = mannwhitneyu(control_guesses, guesses, alternative='two-sided')
                if 0.01 < pvalue <= 0.05:
                    lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(guesses))+')*')
                elif 0.001 < pvalue <= 0.01:
                    lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(guesses))+')**')
                elif pvalue <= 0.001:
                    lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(guesses))+')***')
                else:
                    lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(guesses))+')')

                # calculate quartiles and deciles
                quan_min = np.log(np.quantile(guesses, 0.25))
                quan_max = np.log(np.quantile(guesses, 0.75))
                dec_min = np.log(np.quantile(guesses, 0.1))
                dec_max = np.log(np.quantile(guesses, 0.9))

                # plot the thread stats for the given number of views
                axarr[true_dots_idx,idm].plot(np.log(np.median(guesses)),
                                     [view_idx], 'o', markersize=9, c=colors[true_dots_idx])
                axarr[true_dots_idx,idm].plot(np.log(np.mean(guesses)),
                                     [view_idx], 'v', markersize=7, c=colors[true_dots_idx])
                axarr[true_dots_idx,idm].hlines(y=view_idx, xmin=quan_min, xmax=quan_max,
                                       linewidth=4, color=colors[true_dots_idx])
                axarr[true_dots_idx,idm].hlines(y=view_idx, xmin=dec_min, xmax=dec_max,
                                       linewidth=2, color=colors[true_dots_idx])

            axarr[true_dots_idx,idm].axvline(x=np.log(true_number_of_dots),
                                   linewidth=2, color='#34495e')


            # plotting paraphernalia for the subfigure for the number of dots
            if idm == 0:
                axarr[true_dots_idx,idm].set_yticks([-0.5, 0, 1, 2, 3, 3.5])
                axarr[true_dots_idx,idm].tick_params(axis="y", length=0)  # remove small ticks
                axarr[true_dots_idx,idm].set_yticklabels(lv[true_dots_idx], fontsize=14, ha='left')
                axarr[true_dots_idx,idm].tick_params(axis='y', direction='out', pad=120)

                # make secondary yaxis in order to write image label
                axr = axarr[true_dots_idx,idm].twinx()
                axr.set_yticks([])
                axr.set_ylabel('d='+str(true_number_of_dots), fontsize=18)
                # axr.tick_params(axis='y', direction='out',pad=120)
            else:
                # make secondary yaxis in order to write image label
                axrr = axarr[true_dots_idx,idm].twinx()
                axrr.set_yticks([-0.5, 0, 1, 2, 3, 3.5])
                axrr.set_yticklabels(lv[true_dots_idx], fontsize=14, ha='left')
                axrr.tick_params(axis='y', length=0)
                # axrr.tick_params(axis='y', direction='out', pad=120)

                axarr[true_dots_idx,idm].set_yticks([-0.5, 0, 1, 2, 3, 3.5])
                # axarr.tick_params(axis="y", length=0)  # remove small ticks
                axarr[true_dots_idx,idm].set_yticklabels([])
            # plt.yticks([])  # removing axis numbers with empty list

        for i in range(4):
            axarr[i,idm].tick_params(axis="both", length=0)
            axarr[i,idm].grid(True)

        axarr[3,idm].set_xticks([3.5,4,5,6,7,8,9])
        axarr[3,idm].set_xticklabels(['','4.0','5.0','6.0','7.0','8.0','9.0'], fontdict={'fontsize': 14})
        axarr[3,idm].set_xlim([3.7, 9.8])  # set x-axis dimensions

    axarr[3,0].set_xlabel('log(estimate)', fontsize=18)
    axarr[3,1].set_xlabel('log(estimate)', fontsize=18)
    axarr[0][0].set_title('Historical threads', fontsize='xx-large')
    axarr[0][1].set_title('Manipulated threads', fontsize='xx-large')
    # plt.subplots_adjust(wspace=0.1)
    plt.tight_layout()

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    plt.savefig(os.path.join(PLOTS_DIR, 'summary_stats_plot.png'), transparent=True, dpi=300)
    plt.savefig(os.path.join(PLOTS_DIR, 'summary_stats_plot.pdf'), transparent=True, dpi=300)
    plt.show()


# load the data
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_excel(datafile))
    df_all = df_all.append(df)
#df_all = pd.DataFrame(pd.read_excel(datafiles))
plot_stats(df_all)
