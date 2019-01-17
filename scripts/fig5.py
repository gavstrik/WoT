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
# plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"

"""
plot the social influence of influencers and non-influencers.
Use the untrimmed data and trim after having calculated the social
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

    ids = df['id'].values
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
            if not g in si:
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
    data, _ = zip(*influence_tup)
    max_error_rate = 10
    # lower and upper bounds
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    outliers = [i for i in influence_tup if i[0] >= lower and i[0] <= upper]
    # remove outliers
    data = [i for i in influence_tup if i[0] >= lower and i[0] <= upper]
    return data


def plot_MRE(df_all):
    sessions = ['3cda8qpn', 'bvpj37io', 'ck291lk5', '2hxe3g0w', '5du4txa7',
                'spw8qdcd', 'hhb0if6e', 'wv4xujg7', 'nuf4ogmf', 'i0x2cimn']
    colors = ['#d94801', '#6a51a3']

    MREs = []
    MREs_i_CI = []
    MREs_o_CI = []
    for session in sessions:
        df = df_all[df_all.session == session]
        true_number_of_dots = df['dots'].unique().item()
        views = df['views'].unique().item()
        guesses = (np.array(df['guess'].values))
        influence = find_social_influence(df)

        # make tuples of guesses and their influences
        inf_tup = [(i, influence[pos]) for pos, i in enumerate(guesses)]

        # remove outliers (only after influence score is calculated),
        inf_tup = remove_outliers_error_rate(inf_tup, true_number_of_dots)

        # sort, split
        tup_sorted = sorted(inf_tup, key=lambda tup: tup[1])
        low, high = split_list(tup_sorted)
        guess_influencers, _ = zip(*high)
        guess_others, _ = zip(*low)

        # find mean relative errors and append:
        MRE_i = MRE(true_number_of_dots, np.array(guess_influencers))
        MRE_o = MRE(true_number_of_dots, np.array(guess_others))
        MREs.append((MRE_i, MRE_o))

        # find confidence intervals and append
        errors_influencers = RE(true_number_of_dots, np.array(guess_influencers))
        errors_others = RE(true_number_of_dots, np.array(guess_others))
        CI_influencers = sms.DescrStatsW(errors_influencers).tconfint_mean()
        CI_others = sms.DescrStatsW(errors_others).tconfint_mean()
        MREs_i_CI.append(CI_influencers)
        MREs_o_CI.append(CI_others)

    # unpack and define confidence intervals
    error_inf, error_ninf = zip(*MREs)
    inf_lower, inf_upper = zip(*MREs_i_CI)
    inf_lower = np.array(error_inf) - np.array(inf_lower)
    inf_upper = np.array(inf_upper) - np.array(error_inf)
    ninf_lower, ninf_upper = zip(*MREs_o_CI)
    ninf_lower = np.array(error_ninf) - np.array(ninf_lower)
    ninf_upper = np.array(ninf_upper) - np.array(error_ninf)

    # plot
    fig, ax = plt.subplots(1,1)
    x = np.arange(float(len(sessions)))
    ax.errorbar(x, error_inf, yerr=[inf_lower,inf_upper], capsize=5, fmt='o', c=colors[1], label='high-influencers')
    ax.errorbar(x, error_ninf, yerr=[ninf_lower,ninf_upper], capsize=5, fmt='o',c=colors[0], label='low-influencers')

    # plotting paraphernalia
    x_axis_labels = ['55 dots\n3 views', '55 dots\n9 views',
                     '148 dots\n3 views', '148 dots\n9 views',
                     '403 dots\n3 views', '403 dots\n9 views',
                     '1097 dots\n3 views', '1097 dots\n9 views',
                     '1233 kilo\n3 views', '1233 kilo\n9 views',
                     '403 dot\n3 views\n(max)']
    plt.xticks([i for i in range(len(sessions))], x_axis_labels, rotation=90)
    ax.yaxis.set_major_formatter(FuncFormatter('{0:.0%}'.format))
    plt.ylabel('Mean relative percentage error')
    legend = ax.legend(loc='upper left') #title='Mean estimation error', loc='upper left')
    plt.tight_layout()
    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig('../plots/Fig5.png', transparent=True, dpi=300)
    plt.savefig('../plots/Fig5.pdf', transparent=True, dpi=300)
    plt.show()
    print('average difference btw high-influencers and low-influencers =',
          np.mean(np.array(error_ninf)-np.array(error_inf)))


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)

plot_MRE(df_all)
