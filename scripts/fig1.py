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
            '../data/dots/all_dots_untrimmed_anonymous.csv',
            '../data/ox/all_ox_untrimmed_anonymous.csv',
            ]


# outliers with an error rate above 10 are removed
def remove_outliers(data, true_number_of_dots):
    max_error_rate = 10
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    newdata = [i for i in data if i >= lower and i <= upper]
    return np.array(newdata)


def plot_stats(df_all):
    dots = [1233, 1097, 403, 148, 55]
    views = [0, 1, 3, 9]
    colors = ['#fdae6b', '#f16913', '#d94801', '#7f2704']

    # global figure settings:
    fig, axarr = plt.subplots(nrows=5, sharex=True, figsize=(11.4, 7))
    plt.xlabel('log(estimate)', fontsize=18)
    # plt.rcParams.update({"yticklabel.direction" : "out"})

    # plot
    lv = [[''] for i in range(5)]
    for true_dots_idx, true_dots in enumerate(dots):
        df_i = df_all[df_all['dots'] == true_dots]
        df_control = df_i[df_i['views'] == 0]
        control_guesses = remove_outliers(df_control.guess.values, true_dots)
        true_number_of_dots = true_dots
        for view_idx, view in enumerate(views):
            df = df_i[df_i['views'] == view]
            guesses = remove_outliers(df.guess.values, true_dots)

            # calculate the p-value using the two-sided Mann-Whitney U test
            _, pvalue = mannwhitneyu(control_guesses, guesses)
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
            axarr[true_dots_idx].plot(np.log(np.median(guesses)),
                                 [view_idx], 'o', markersize=9, c=colors[view_idx])
            axarr[true_dots_idx].plot(np.log(np.mean(guesses)),
                                 [view_idx], 'v', markersize=7, c=colors[view_idx])
            axarr[true_dots_idx].hlines(y=view_idx, xmin=quan_min, xmax=quan_max,
                                   linewidth=4, color=colors[view_idx])
            axarr[true_dots_idx].hlines(y=view_idx, xmin=dec_min, xmax=dec_max,
                                   linewidth=2, color=colors[view_idx])

        axarr[true_dots_idx].axvline(x=np.log(true_number_of_dots),
                               linewidth=2, color='#34495e')

        # plotting paraphernalia for the subfigure for the number of dots
        axarr[true_dots_idx].set_yticks([-0.5, 0, 1, 2, 3, 3.5])
        axarr[true_dots_idx].tick_params(axis="y", length=0)  # remove small ticks
        axarr[true_dots_idx].set_yticklabels(lv[true_dots_idx], fontsize=14, ha='left')
        axarr[true_dots_idx].tick_params(axis='y', direction='out', pad=120)

        # make secondary yaxis in order to write image label
        axr = axarr[true_dots_idx].twinx()
        axr.set_ylabel('d='+str(true_number_of_dots), fontsize=18)
        plt.yticks([])  # removing axis numbers with empty list

    for i in range(5):
        axarr[i].tick_params(axis="x", length=0)
        axarr[i].grid(True)

    axarr[4].set_xticklabels(['','4.0','4.5','5.0','5.5','6.0','6.5','7.0','7.5',''], fontdict={'fontsize': 14})


    axarr[0].set_xlim([3.7, 7.8])  # set x-axis dimensions
    plt.tight_layout()

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    plt.savefig(os.path.join(PLOTS_DIR, 'fig1.png'), transparent=True, dpi=300)
    plt.savefig(os.path.join(PLOTS_DIR, 'fig1.pdf'), transparent=True, dpi=300)
    plt.show()


# load the data
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)

# only plot the data for the history sessions
df_all = df_all[df_all.method == 'history']
plot_stats(df_all)
