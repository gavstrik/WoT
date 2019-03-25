import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.api as sms
import statsmodels.stats.proportion as smsp

plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'

"""
Plotting aggregate measures showing that individuals as well as
the collective becomes more accurate with increasing social information.
"""

datafiles = [
            '../data/dots.csv',
            '../data/ox.csv',
            ]


# outliers with an error rate above 10 are removed
def remove_outliers(data, true_number_of_dots):
    max_error_rate = 10
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    newdata = [i for i in data if i >= lower and i <= upper]
    return np.array(newdata)


# relative error
def RE(aggregate_prediction, truth):
    return abs(truth - aggregate_prediction)/truth


# mean relative error
def MRE(individual_predictions, truth):
    return np.mean(abs(truth - individual_predictions))/truth


# median relative error
def MdRE(individual_predictions, truth):
    return np.median(abs(truth - individual_predictions))/truth


# confidence interval of a median
def medianCI(data, ci, percentile):
    """ This piece of code is stolen from
    https://github.com/minddrummer/median-confidence-interval/blob/master/Median_CI.py#L6
    """
    if type(data) is pd.Series or type(data) is pd.DataFrame:
        data = data.values
    data = np.sort(data)
    N = data.shape[0]
    lowCount, upCount = stats.binom.interval(ci, N, percentile, loc=0)
	#given this: https://onlinecourses.science.psu.edu/stat414/node/316
    #lowCount and upCount both refers to  W's value, W follows binomial Dis.
    #lowCount need to change to lowCount-1, upCount no need to change in python indexing
    lowCount -= 1
    return data[int(lowCount)], data[int(upCount)]


def plot_aggregates(df_alld):
    fig, axes = plt.subplots(nrows=2, ncols=4,
                             sharey=True,
                             figsize=(11.4,9))
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
    trimmed_observations_fixed_v = [0,0,0,0]
    BigN = 0
    for dots_idx, dots in enumerate([55, 148, 403, 1097, 1233]):
        col_error_median = []
        col_error_median_CIs = []
        col_error_mean = []
        col_error_mean_CIs = []
        median_rel_error = []
        median_rel_error_CIs = []
        mean_rel_error = []
        mean_rel_error_CIs = []
        better_than_median = []
        better_than_median_CIs = []
        better_than_mean = []
        better_than_mean_CIs = []
        bonus = []
        bonus_CIs = []
        df_all = df_alld[df_alld.d == dots]
        N = 0
        for position, views in enumerate([0, 1, 3, 9]):
            df = df_all[df_all.v == views]
            untrimmed_guesses = df.guess.values
            guesses = remove_outliers(untrimmed_guesses, dots)
            N += len(guesses)
            trimmed_observations_fixed_v[position] += len(guesses)

            # find the collective error of the median:
            error_of_median = RE(np.median(guesses), dots)
            col_error_median.append(error_of_median)
            avg_col_error_median[position] += error_of_median
            # find error of the confidence interval of the median
            error_of_lower_ci_median = RE(medianCI(guesses, 0.95, 0.50)[0], dots)
            error_of_upper_ci_median = RE(medianCI(guesses, 0.95, 0.50)[1], dots)
            col_error_median_CIs.append((error_of_lower_ci_median,error_of_upper_ci_median))

            # find the collective error of the mean:
            error_of_mean = RE(np.mean(guesses), dots)
            col_error_mean.append(error_of_mean)
            avg_col_error_mean[position] += error_of_mean
            # find the error of the confidence interval of the mean
            error_of_lower_ci_mean = RE(sms.DescrStatsW(guesses).tconfint_mean()[0], dots)
            error_of_upper_ci_mean = RE(sms.DescrStatsW(guesses).tconfint_mean()[1], dots)
            col_error_mean_CIs.append((error_of_lower_ci_mean, error_of_upper_ci_mean))

            # find the median individual error:
            median_rel_error.append(MdRE(guesses, dots))
            avg_median_rel_error[position] += MdRE(guesses, dots)
            # find the confidence interval of the median relative error:
            relative_errors = RE(guesses, dots)
            ci_of_median_rel_error = medianCI(relative_errors, 0.95, 0.50)
            median_rel_error_CIs.append(ci_of_median_rel_error)

            # find the mean individual error:
            mean_rel_error.append(MRE(guesses, dots))
            avg_mean_rel_error[position] += MRE(guesses, dots)
            # find the confidence interval of the mean relative error:
            relative_errors = RE(guesses, dots)
            ci_of_mean_rel_error = sms.DescrStatsW(relative_errors).tconfint_mean()
            mean_rel_error_CIs.append(ci_of_mean_rel_error)

            # find the fraction of individuals better than the median and mean:
            n_median = [1 for g in guesses
                        if RE(g, dots) < error_of_median]
            n_mean = [1 for g in guesses
                      if RE(g, dots) < error_of_mean]
            better_than_median.append(len(n_median)/len(guesses))
            better_than_mean.append(len(n_mean)/len(guesses))
            avg_better_than_median[position] += (len(n_median)/len(guesses))
            avg_better_than_mean[position] += (len(n_mean)/len(guesses))
            # find the confidence intervals
            better_than_median_CIs.append(smsp.proportion_confint(len(n_median), len(guesses), alpha=0.05, method='normal'))
            better_than_mean_CIs.append(smsp.proportion_confint(len(n_mean), len(guesses), alpha=0.05, method='normal'))

            # find the fraction of participants winning the bonus:
            bonus.append(df['bonus'].sum()/len(guesses))
            avg_bonus[position] += df['bonus'].sum()/len(guesses)
            bonus_CIs.append(smsp.proportion_confint(df['bonus'].sum(), len(guesses), alpha=0.05, method='normal'))

        x = np.arange(4) # scalar for x-axis in errorbars
        # plot the collective error of the median and associated CIs:
        col_error_median_lower, col_error_median_upper = zip(*col_error_median_CIs)
        col_error_median_lower = np.array(col_error_median) - np.array(col_error_median_lower)
        col_error_median_upper = np.array(col_error_median_upper) - np.array(col_error_median)
        axes[0].errorbar(x, col_error_median, yerr=[col_error_median_lower, col_error_median_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8, capsize=3, linestyle='-',
                       linewidth=2, fmt='o', c=colors[dots_idx], label='d='+str(dots))

        # plot the median relative error and associated CIs:
        median_rel_error_lower, median_rel_error_upper = zip(*median_rel_error_CIs)
        median_rel_error_lower = np.array(median_rel_error) - np.array(median_rel_error_lower)
        median_rel_error_upper = np.array(median_rel_error_upper) - np.array(median_rel_error)
        axes[1].errorbar(x, median_rel_error, yerr=[median_rel_error_lower, median_rel_error_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8, capsize=3, linestyle='-',
                       linewidth=2, fmt='o', c=colors[dots_idx], label='d='+str(dots))

        # plot the fraction of participants better than the median and the CIs:
        better_than_median_lower, better_than_median_upper = zip(*better_than_median_CIs)
        better_than_median_lower = np.array(better_than_median) - np.array(better_than_median_lower)
        better_than_median_upper = np.array(better_than_median_upper) - np.array(better_than_median)
        axes[2].errorbar(x, better_than_median, yerr=[better_than_median_lower, better_than_median_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8, capsize=3, linestyle='-',
                       linewidth=2, fmt='o', c=colors[dots_idx], label='d='+str(dots))

        # plot the fraction getting a bonus and associated CIs:
        bonus_lower, bonus_upper = zip(*bonus_CIs)
        bonus_lower = np.array(bonus) - np.array(bonus_lower)
        bonus_upper = np.array(bonus_upper) - np.array(bonus)
        axes[3].errorbar(x, bonus, yerr=[bonus_lower, bonus_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8, capsize=3, linestyle='-',
                       linewidth=2, fmt='o', c=colors[dots_idx], label='d='+str(dots)+' (N='+str(N)+')')

        # plot the collective error of the mean and associated CIs:
        col_error_mean_lower, col_error_mean_upper = zip(*col_error_mean_CIs)
        col_error_mean_lower = np.array(col_error_mean) - np.array(col_error_mean_lower)
        col_error_mean_upper = np.array(col_error_mean_upper) - np.array(col_error_mean)
        axes[4].errorbar(x, col_error_mean, yerr=[col_error_mean_lower, col_error_mean_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8, capsize=3, linestyle='-',
                       linewidth=2, fmt='o', c=colors[dots_idx], label='d='+str(dots))

        # plot the mean relative error and associated CIs:
        mean_rel_error_lower, mean_rel_error_upper = zip(*mean_rel_error_CIs)
        mean_rel_error_lower = np.array(mean_rel_error) - np.array(mean_rel_error_lower)
        mean_rel_error_upper = np.array(mean_rel_error_upper) - np.array(mean_rel_error)
        axes[5].errorbar(x, mean_rel_error, yerr=[mean_rel_error_lower, mean_rel_error_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8, capsize=3, linestyle='-',
                       linewidth=2, fmt='o', c=colors[dots_idx], label='d='+str(dots))

        # plot the fraction of participants better than the mean and the CIs:
        better_than_mean_lower, better_than_mean_upper = zip(*better_than_mean_CIs)
        better_than_mean_lower = np.array(better_than_mean) - np.array(better_than_mean_lower)
        better_than_mean_upper = np.array(better_than_mean_upper) - np.array(better_than_mean)
        axes[6].errorbar(x, better_than_mean, yerr=[better_than_mean_lower, better_than_mean_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8, capsize=3, linestyle='-',
                       linewidth=2, fmt='o', c=colors[dots_idx], label='d='+str(dots))

    BigN = np.sum(trimmed_observations_fixed_v)
    # plot the averages across all d:
    axes[0].plot(np.array(avg_col_error_median)/5, zorder=6, marker='o', linestyle='-', linewidth=3, markersize=8, color=colors[dots_idx+1], label='average')
    axes[1].plot(np.array(avg_median_rel_error)/5, zorder=6, marker='o', linestyle='-', linewidth=3, markersize=8, color=colors[dots_idx+1], label='average')
    axes[2].plot(np.array(avg_better_than_median)/5, zorder=6, marker='o', linestyle='-', linewidth=3, markersize=8, color=colors[dots_idx+1], label='average')
    axes[3].plot(np.array(avg_bonus)/5, zorder=6, marker='o', linestyle='-', linewidth=3, markersize=8, color=colors[dots_idx+1], label='average')
    axes[4].plot(np.array(avg_col_error_mean)/5, zorder=6, marker='o', linestyle='-', linewidth=3, markersize=8, color=colors[dots_idx+1], label='average')
    axes[5].plot(np.array(avg_mean_rel_error)/5, zorder=6, marker='o', linestyle='-', linewidth=3, markersize=8, color=colors[dots_idx+1], label='average')
    axes[6].plot(np.array(avg_better_than_mean)/5, zorder=6, marker='o', linestyle='-', linewidth=3, markersize=8, color=colors[dots_idx+1], label='average')
    print('average bonuses:', np.array(avg_bonus)/5)

    # plotting paraphernalia
    handles, labels = axes[3].get_legend_handles_labels()
    leg = plt.figlegend(handles, labels, loc=(0.77,0.3), prop={'size': 13})
    # leg.set_title("Legend", prop = {'size':'x-large'})
    txtA = fig.text(0.162, .99, 'A', fontsize='xx-large', fontweight='bold', ha='center')
    txtB = fig.text(0.402, .99, 'B', fontsize='xx-large', fontweight='bold', ha='center')
    txtC = fig.text(0.642, .99, 'C', fontsize='xx-large', fontweight='bold', ha='center')
    txtD = fig.text(0.885, .99, 'D', fontsize='xx-large', fontweight='bold', ha='center')
    txtv = fig.text(0.525, .0, 'v', fontsize='xx-large', fontweight='bold', ha='center')

    # set y-labels
    axes[0].set_ylabel('Collective error of median', fontsize='xx-large')
    axes[1].set_ylabel('Median individual error', fontsize='xx-large')
    axes[2].set_ylabel('Better than median', fontsize='xx-large')
    axes[3].set_ylabel('Fraction winning a bonus', fontsize='xx-large')
    axes[4].set_ylabel('Collective error of mean', fontsize='xx-large')
    axes[5].set_ylabel('Mean individual error', fontsize='xx-large')
    axes[6].set_ylabel('Better than mean', fontsize='xx-large')
    # remove small ticks
    for i in [1,2,3,5,6]:
        axes[i].tick_params(axis="y", length=0)
    for i in range(3):
        axes[i].set_xticklabels('')
        axes[i].tick_params(axis="x", length=0)
    # rename xticks
    x_axis_labels = ['', '0', '1', '3', '9', '']
    for i in range(3, 7):
        axes[i].set_xticklabels(x_axis_labels, fontdict={'fontsize': 14})
    # resize yticks
    axes[0].set_yticklabels(['','0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'], fontdict={'fontsize': 14})
    axes[4].set_yticklabels(['','0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'], fontdict={'fontsize': 14})

    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    fig.savefig(os.path.join(PLOTS_DIR, 'fig2.png'), transparent=True,
                bbox_extra_artists=(txtA,txtB,txtC,txtD,txtv),
                bbox_inches='tight', dpi=300)
    fig.savefig(os.path.join(PLOTS_DIR, 'fig2.pdf'), transparent=True,
                bbox_extra_artists=(txtA,txtB,txtC,txtD,txtv),
                bbox_inches='tight', dpi=300)
    plt.show()


# load the data
dataframe = pd.DataFrame()
for datafile in datafiles:
    dataframe = dataframe.append(pd.DataFrame(pd.read_csv(datafile)), sort=True)
dataframe = dataframe[dataframe.method == 'history']

plot_aggregates(dataframe)
