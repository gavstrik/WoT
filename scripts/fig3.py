import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"

PLOTS_DIR = '../plots'

"""
Plotting aggregate measures showing that individuals as well as
the collective becomes more accurate with increasing social information.
"""

datafiles = [
            '../data/dots/all_dots_trimmed_anonymous.csv',
            '../data/ox/all_ox_untrimmed_anonymous.csv',
            ]


# relative error
def RE(aggregate_prediction, truth):
    return abs(truth - aggregate_prediction)/truth


# mean relative error
def MRE(individual_predictions, truth):
    return np.mean(abs(truth - individual_predictions))/truth


# median relative error
def MdRE(individual_predictions, truth):
    return np.median(abs(truth - individual_predictions))/truth


def plot_aggregates(df_all):
    fig, axes = plt.subplots(nrows=2, ncols=4,
                             sharex=True, sharey=True,
                             figsize=(6,6))
    axes[-1, -1].axis('off')  # do not show the last subplot
    axes = axes.flatten()  # needs to be flattened for iteration
    colors = ['#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#3f007d', '#d94801']

    # plot aggregate measures as a function of v:
    avg_col_error_median = [0,0,0,0]
    avg_col_error_mean = [0,0,0,0]
    avg_median_rel_error = [0,0,0,0]
    avg_mean_rel_error = [0,0,0,0]
    avg_better_than_median = [0,0,0,0]
    avg_better_than_mean = [0,0,0,0]
    avg_bonus = [0,0,0,0]

    # NB: 1233 is the ox experiment
    for dots_idx, dots in enumerate([55, 148, 403, 1097, 1233]):
        col_error_median = []
        col_error_mean = []
        median_rel_error = []
        mean_rel_error = []
        better_than_median = []
        better_than_mean = []
        bonus = []
        for position, views in enumerate([0, 1, 3, 9]):
            df = df_all[(df_all['dots'] == dots) &
                        (df_all['views'] == views)]

            # Skip the treatment combinations that don't exist
            # (e.g. 'min'-method for dots experiments)
            if len(df) == 0:
                continue

            # find the collective error of the median:
            median_error = RE(df.guess.median(), df.dots.unique().item())
            col_error_median.append(median_error)
            avg_col_error_median[position] += median_error/5

            # find the collective error of the mean:
            mean_error = RE(df.guess.mean(), df.dots.unique().item())
            col_error_mean.append(mean_error)
            avg_col_error_mean[position] += mean_error/5

            # find the median individual error:
            median_rel_error.append(MdRE(df.guess.values,
                                         df.dots.unique().item()))
            avg_median_rel_error[position] += MdRE(df.guess.values,
                                           df.dots.unique().item())/5

            # find the mean individual error:
            mean_rel_error.append(MRE(df.guess.values, df.dots.unique().item()))
            avg_mean_rel_error[position] += MRE(df.guess.values, df.dots.unique().item())/5

            # find the fraction of individuals better than the median and mean:
            n_median = [1 for g in df.guess.values
                        if RE(g, df.dots.unique().item()) < median_error]
            n_mean = [1 for g in df.guess.values
                      if RE(g, df.dots.unique().item()) < mean_error]
            better_than_median.append(len(n_median)/len(df.guess.values))
            better_than_mean.append(len(n_mean)/len(df.guess.values))
            avg_better_than_median[position] += (len(n_median)/len(df.guess.values))/5
            avg_better_than_mean[position] += (len(n_mean)/len(df.guess.values))/5

            # find the fraction of participants winning the bonus:
            bonus.append(df['bonus'].sum()/len(df.guess.values))
            avg_bonus[position] += df['bonus'].sum()/len(df.guess.values)/5

        # plot the upper/lower row of measures for the median/mean:
        axes[0].plot(col_error_median, marker='o', linestyle='-', color=colors[dots_idx], label='d='+str(dots))
        axes[1].plot(median_rel_error, marker='o', linestyle='-', color=colors[dots_idx], label='d='+str(dots))
        axes[2].plot(better_than_median, marker='o', linestyle='-', color=colors[dots_idx], label='d='+str(dots))
        axes[3].plot(bonus, marker='o', linestyle='-', color=colors[dots_idx], label='d='+str(dots))
        axes[4].plot(col_error_mean, marker='o', linestyle='-', color=colors[dots_idx], label='d='+str(dots))
        axes[5].plot(mean_rel_error, marker='o', linestyle='-', color=colors[dots_idx], label='d='+str(dots))
        axes[6].plot(better_than_mean, marker='o', linestyle='-', color=colors[dots_idx], label='d='+str(dots))

    # plot the averages across all d:
    axes[0].plot(avg_col_error_median, marker='o', linestyle='-', linewidth=3, color=colors[dots_idx+1], label='average')
    axes[1].plot(avg_median_rel_error, marker='o', linestyle='-', linewidth=3, color=colors[dots_idx+1], label='average')
    axes[2].plot(avg_better_than_median, marker='o', linestyle='-', linewidth=3, color=colors[dots_idx+1], label='average')
    axes[3].plot(avg_bonus, marker='o', linestyle='-', linewidth=3, color=colors[dots_idx+1], label='average')
    axes[4].plot(avg_col_error_mean, marker='o', linestyle='-', linewidth=3, color=colors[dots_idx+1], label='average')
    axes[5].plot(avg_mean_rel_error, marker='o', linestyle='-', linewidth=3, color=colors[dots_idx+1], label='average')
    axes[6].plot(avg_better_than_mean, marker='o', linestyle='-', linewidth=3, color=colors[dots_idx+1], label='average')

    # plotting paraphernalia
    handles, labels = axes[3].get_legend_handles_labels()
    plt.figlegend(handles, labels, loc=(0.795,0.25), title='legend', ncol=1)
    txtA = fig.text(0.195, .99, 'A', fontsize='x-large', ha='center')
    txtB = fig.text(0.425, .99, 'B', fontsize='x-large', ha='center')
    txtC = fig.text(0.655, .99, 'C', fontsize='x-large', ha='center')
    txtD = fig.text(0.885, .99, 'D', fontsize='x-large', ha='center')
    txtv = fig.text(0.54, .0, 'v', fontsize='large', ha='center')

    x_axis_labels = ['0', '1', '3', '9']  # rename xticks
    plt.xticks([i for i in range(4)], x_axis_labels)

    # set y-labels
    axes[0].set_ylabel('Collective error of median')
    axes[1].set_ylabel('Median individual error')
    axes[2].set_ylabel('Better than median')
    axes[3].set_ylabel('Fraction winning a bonus')
    axes[4].set_ylabel('Collective error of mean')
    axes[5].set_ylabel('Mean individual error')
    axes[6].set_ylabel('Better than mean')

    # remove small ticks
    for i in range(7):
        axes[i].tick_params(axis="both", length=0)
        axes[i].grid(True)
    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    output_path = os.path.join(PLOTS_DIR, 'fig3.png')

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    fig.savefig(output_path, transparent=True, bbox_extra_artists=(txtA,txtB,txtC,txtD,txtv),
                bbox_inches='tight', dpi=300)
    fig.savefig(output_path, transparent=True, bbox_extra_artists=(txtA,txtB,txtC,txtD,txtv),
                bbox_inches='tight', dpi=300)
    plt.show()


# load the data
dataframe = pd.DataFrame()
for datafile in datafiles:
    dataframe = dataframe.append(pd.DataFrame(pd.read_csv(datafile)))

# only plot the data for the history sessions
dataframe = dataframe[(dataframe['method'] == 'history')]
plot_aggregates(dataframe)
