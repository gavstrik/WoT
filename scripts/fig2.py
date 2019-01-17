import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, gmean
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"

"""
This python script plotting summary statistics for
the dots-experiments on Amazon Mechanical Turk.
"""

data_file = pd.read_csv('../data/dots/all_dots_trimmed_anonymous.csv')


def plot_stats(df):
    dots = [55, 148, 403, 1097]
    views = [0, 1, 3, 9]
    colors = ['#fdae6b', '#f16913', '#d94801', '#7f2704']

    # # global figure settings:
    fig, axarr = plt.subplots(nrows=4, sharex=True, figsize=(5.5, 4))
    plt.xlabel('log(estimate)', fontweight='bold')

    # plot
    lv = [[''] for i in range(5)]
    for position, true_number_of_dots in enumerate(dots):
        df_i = df_all[df_all['dots'] == true_number_of_dots]
        for pos, view in enumerate(views):
            df = df_i[df_i['views'] == view]
            sessions = df.session.unique()

            if view == 0:
                control = df.guess.values

            # calculate the p-value using the two-sided Mann-Whitney U test
            _, pvalue = mannwhitneyu(control, df.guess.values)
            if pvalue <= 0.05 and pvalue > 0.01:
                lv[position].append('v='+str(view)+' (N='+str(len(df))+')*   ')
            elif pvalue <= 0.01 and pvalue > 0.001:
                lv[position].append('v='+str(view)+' (N='+str(len(df))+')** ')
            elif pvalue <= 0.001:
                lv[position].append('v='+str(view)+' (N='+str(len(df))+')***')
            else:
                lv[position].append('v='+str(view)+' (N='+str(len(df))+')    ')

            # calculate quartiles and deciles
            quan_min = df['guess'].apply(np.log).quantile(0.25)
            quan_max = df['guess'].apply(np.log).quantile(0.75)
            dec_min = df['guess'].apply(np.log).quantile(0.1)
            dec_max = df['guess'].apply(np.log).quantile(0.9)

            log_guesses = df['guess'].apply(np.log).values

            # plot
            axarr[position].plot(np.log(np.median(df['guess'].values)),
                                 [pos], 'o', markersize=9, c=colors[pos])
            axarr[position].plot(np.log(np.mean(df['guess'].values)),
                                 [pos], 'v', markersize=7, c=colors[pos])
            axarr[position].hlines(y=pos, xmin=quan_min, xmax=quan_max,
                                   linewidth=4, color=colors[pos])
            axarr[position].hlines(y=pos, xmin=dec_min, xmax=dec_max,
                                   linewidth=2, color=colors[pos])
        axarr[position].axvline(x=np.log(true_number_of_dots),
                               linewidth=1, color='black')

        # plotting paraphernalia
        axarr[position].set_yticks([-0.5,0,1,2,3,3.5])
        axarr[position].tick_params(axis="y", length=0)  # remove small ticks
        axarr[position].set_yticklabels(lv[position])
        axr = axarr[position].twinx()
        axr.set_ylabel('d='+str(true_number_of_dots), fontweight='bold')
        plt.yticks([])  # removing axis numbers with empty list

    for i in range(4):
        axarr[i].tick_params(axis="x", length=0)
        axarr[i].grid(True)

    axarr[0].set_xlim([3.7, 7.7])  # set x-axis dimensions
    plt.tight_layout()
    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig('../plots/Fig2.png', transparent=True, dpi=300)
    plt.savefig('../plots/Fig2.pdf', transparent=True, dpi=300)
    plt.show()


# main code
df_all = pd.DataFrame(data_file)

# use only history data
df_all = df_all[df_all.method == 'history']
plot_stats(df_all)
