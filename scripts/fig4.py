import os
import numpy as np
import pandas as pd
import statistics
import statsmodels.stats.api as sms
from numpy import percentile
import scipy.stats as st
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import style
from matplotlib.ticker import FuncFormatter
import matplotlib
plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'

"""
plot the social influence of influencers and non-influencers.
Use the raw data and trim after having calculated the social
influence score for each participant.
"""

datafiles = [
            '../data/dots/all_dots_untrimmed_anonymous.csv',
            '../data/ox/all_ox_untrimmed_anonymous.csv',
            ]
noise = 0.01


def social_influence_scores(g, guess_vector):
    influences = 1/np.array([abs(g-h) + noise for h in guess_vector])
    normalizer = sum(influences)
    scores = influences/normalizer
    return scores


# function for finding nearest guess
def find_nearest(array, value):
    array = np.asarray(array)
    # because argmin() always returns the first number in the case of multiple
    # minimum values, we are good with the code as is:
    idx = (np.abs(array - value)).argmin()
    return (idx, array[idx])


# relative error
def RE(targets, predictions):
    predictions = np.array(predictions)
    return abs((targets - predictions)/targets)


# mean relative error
def MRE(targets, predictions):
    predictions = np.array(predictions)
    return np.mean(abs((targets - predictions)/targets))


def find_social_influence(df):
    # initialize dict for social influence
    si = {}

    ids = df.id.values

    for id in ids:
        own_guess = df[df['id'] == id]['guess'].item()
        seen_ids = df[df['id'] == id]['hist'].item()
        seen_ids = pd.eval(seen_ids)
        seen_ids = [h for h in seen_ids if h in ids]
        if not seen_ids:
            continue
        seen_guesses = [df[df['id'] == g]['guess'].item() for g in seen_ids]

        # update the social influence values of the seen ids:
        for pos, g in enumerate(seen_ids):
            if g not in si:
                si[g] = social_influence_scores(own_guess, seen_guesses)[pos]
            else:
                si[g] += social_influence_scores(own_guess, seen_guesses)[pos]

    # make a new list for the social influence of all seen ids
    infl = []
    for id in ids:
        if id in si:
            infl.append(si[id])
        else:
            infl.append(0)

    return infl


def split_list(sorted_list):
    half = int(len(sorted_list)/2)
    return sorted_list[:half], sorted_list[half:]


def remove_outliers_error_rate(influence_tup, true_number_of_dots):
    # data, _ = zip(*influence_tup)
    max_error_rate = 10
    # lower and upper bounds
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    data = [i for i in influence_tup if i[0] >= lower and i[0] <= upper]
    return data


def plot_MRE(df_all):
    colors = ['#d94801', '#6a51a3']
    MREs = []
    MREs_i_CI = []
    MREs_ni_CI = []
    pops = []
    for true_number_of_dots in [55, 148, 403, 1097, 1233]:
        for views in [3, 9]:
            df1 = df_all[(df_all.dots == true_number_of_dots) & (df_all.views == views)]
            sessions = df1.session.unique()
            influence_tuples = []
            for session in sessions:
                df = df1[df1.session == session]
                guesses = (np.array(df.guess.values))
                influence = find_social_influence(df)
                # make a long list of tuples containing guesses and their influences
                influence_tuples.extend([(guess, influence[guess_idx]) for guess_idx, guess in enumerate(guesses)])

            # remove outliers (only after influence score is calculated),
            influence_tuples = remove_outliers_error_rate(influence_tuples, true_number_of_dots)
            pops.append(len(influence_tuples)) # count them

            # sort guesses by ascending social influence
            tup_sorted = sorted(influence_tuples, key=lambda tup: tup[1])
            low, high = split_list(tup_sorted)
            guesses_of_influencers, _ = zip(*high)
            guesses_of_non_influencers, _ = zip(*low)

            # find mean relative errors and append:
            MRE_i = MRE(true_number_of_dots, np.array(guesses_of_influencers))
            MRE_ni = MRE(true_number_of_dots, np.array(guesses_of_non_influencers))
            MREs.append((MRE_i, MRE_ni))

            # find confidence intervals of the mean relative error and append
            errors_influencers = RE(true_number_of_dots, np.array(guesses_of_influencers))
            errors_noninfluencers = RE(true_number_of_dots, np.array(guesses_of_non_influencers))
            CI_influencers = sms.DescrStatsW(errors_influencers).tconfint_mean()
            CI_noninfluencers = sms.DescrStatsW(errors_noninfluencers).tconfint_mean()
            MREs_i_CI.append(CI_influencers)
            MREs_ni_CI.append(CI_noninfluencers)

    # unpack and format confidence intervals
    error_inf, error_ninf = zip(*MREs)
    inf_lower, inf_upper = zip(*MREs_i_CI)
    inf_lower = np.array(error_inf) - np.array(inf_lower)
    inf_upper = np.array(inf_upper) - np.array(error_inf)
    ninf_lower, ninf_upper = zip(*MREs_ni_CI)
    ninf_lower = np.array(error_ninf) - np.array(ninf_lower)
    ninf_upper = np.array(ninf_upper) - np.array(error_ninf)

    # plot
    fig, ax = plt.subplots(1, 1, figsize=(8.7, 6))
    x = np.arange(float(10))
    ax.errorbar(x, error_inf, yerr=[inf_lower, inf_upper], capsize=5, fmt='o',
                elinewidth=2, markeredgewidth=2, markersize=8, c=colors[1], label='high-influencers')
    ax.errorbar(x, error_ninf, yerr=[ninf_lower, ninf_upper], capsize=5,
                elinewidth=2, markeredgewidth=2, markersize=8,
                fmt='o', c=colors[0], label='low-influencers')

    # plotting paraphernalia
    x_axis_labels = [
                     'd=55, v=3\n(N='+str(pops[0])+')',
                     'd=55, v=9\n(N='+str(pops[1])+')',
                     'd=148, v=3\n(N='+str(pops[2])+')',
                     'd=148, v=9\n(N='+str(pops[3])+')',
                     'd=403, v=3\n(N='+str(pops[4])+')',
                     'd=403, v=9\n(N='+str(pops[5])+')',
                     'd=1097, v=3\n(N='+str(pops[6])+')',
                     'd=1097, v=9\n(N='+str(pops[7])+')',
                     'd=1233, v=3\n(N='+str(pops[8])+')',
                     'd=1233, v=9\n(N='+str(pops[9])+')',
                     ]
    plt.xticks([i for i in range(10)], x_axis_labels, fontsize='x-large', rotation=90)
    plt.yticks([0,.2,.4,.6,.8,1], [0,.2,.4,.6,.8,1], fontsize=14)
    plt.ylabel('Mean relative error', fontsize='xx-large')
    ax.legend(loc='upper left', prop={'size': 16})
    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig(os.path.join(PLOTS_DIR, 'fig4.png'), transparent=True, dpi=300)
    plt.savefig(os.path.join(PLOTS_DIR, 'fig4.pdf'), transparent=True, dpi=300)
    plt.show()

    print('average difference btw high-influencers and low-influencers =',
          np.mean(np.array(error_ninf)-np.array(error_inf)))


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)
df_all = df_all[df_all.method == 'history']

plot_MRE(df_all)
