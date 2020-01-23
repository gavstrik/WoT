import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from scipy import stats
import statsmodels.stats.api as sms
from matplotlib.ticker import FuncFormatter
import matplotlib
# plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'

"""
model simulation + calculation of social influence score:
subjects have an initial "hunch" (taken from the control treatment),
see v previous guesses, take the mean of all numbers and report it as
their estimate.
"""

datafiles = [
            '../data/dots.csv',
            '../data/ox.csv',
            ]

# random.seed(4)
tlength = 400
noise = 0.01
simuls = 1


def social_influence_scores(g, guess_vector):
    influences = 1/np.array([abs(g-h) + noise for h in guess_vector])
    normalizer = sum(influences)
    scores = influences/normalizer
    return scores


# relative error
def RE(targets, predictions):
    predictions = np.array(predictions)
    return abs((targets - predictions)/targets)


# mean relative error
def MRE(targets, predictions):
    predictions = np.array(predictions)
    return np.mean(abs((targets - predictions)/targets))


def split_list(sorted_list):
    half = int(len(sorted_list)/2)
    return sorted_list[:half], sorted_list[half:]


# outliers with an error rate above 10 are removed
def remove_outliers(data, true_number_of_dots):
    max_error_rate = 10
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    newdata = [i for i in data if i >= lower and i <= upper]
    return np.array(newdata)


def find_social_influence(thread, histories):
    # initialize dict for social influence
    si = {}

    for position, t in enumerate(thread):
        seen_guesses = histories[position]

        # update the social influence values of the seen guesses:
        for pos, g in enumerate(seen_guesses):
            if not g in si:
                si[g] = social_influence_scores(t, seen_guesses)[pos]
            else:
                si[g] += social_influence_scores(t, seen_guesses)[pos]

    # make a new list for the social influence of all seen ids
    infl = []
    si[thread[-1]] = 0  # the last id has always social influence = 0
    for t in thread:
        infl.append(si[t])

    return infl


def generate_simulation_data(df_all):
    fig, ax = plt.subplots(1,1)
    colors = ['#d94801', '#6a51a3']

    MREs = []
    MREs_i_CI = []
    MREs_o_CI = []
    for position, d in enumerate([55,148,403,1097,1233]):
        # define the control group
        df = df_all[df_all.d == d]
        df_control = df[df.v == 0]
        control = remove_outliers(df_control.guess.values, d)
        # sd = np.std(control)
        sd = d
        normal_dist = list(np.random.normal(loc=d, scale=sd, size=len(control)))
        lognormal_dist = list(np.random.lognormal(np.log(np.mean(control)), 1, size=len(control)))
        lognormal_dist = list(remove_outliers(lognormal_dist, d))

        for v in [3, 9]:
            mi = []
            mo = []
            w_self = 9  # these are some supposed weights by which players may value their own hunch and the estimates by others
            w_others = 1


            # simluate with 'simul' runs:
            for i in range(simuls):
                random_draws = random.sample(list(control), tlength)
                thread = []
                histories = []
                thread.extend(random_draws[:v])  # fill with the first v samples
                histories.extend([random_draws[:v] for i in range(v)])  # fill with the first v histories

                # simulate a thread of length tlength
                for g in range(v, tlength):
                    history = thread[-v:]  # these are the previous estimates
                    own_hunch = random_draws[g]

                    # Here the naive deGroot model of learning:
                    # final_guess = 1 / (v + 1) * (own_hunch + sum(history))

                    # Here some alternative models as suggested by reviewer
                    # final_guess = random.sample(normal_dist, 1)[0] # guess as PNAS reviewer #1 suggested
                    # final_guess = random.sample(lognormal_dist, 1)[0]

                    # Feeding the controls only also gives a sign. diff.:
                    # final_guess = random.sample(list(control), 1)[0]

                    # this variant uses naive deGroot but with weights:
                    # final_guess = (own_hunch * w_self + (sum(history) * w_others))/(w_self + (w_others * v))

                    # This model assumes that players guess randomly around the sample mean:
                    sample_normal_dist = list(np.random.normal(loc=np.mean(history), scale=d/5, size=10))
                    # print(sample_normal_dist)
                    final_guess = random.sample(sample_normal_dist, 1)[0]

                    print(own_hunch, np.mean(history), final_guess)
                    thread.append(final_guess)
                    histories.append(history)

                influence = find_social_influence(thread, histories)

                # make tuples of guesses and their influences
                inf_tup = [(i, influence[pos]) for pos, i in enumerate(thread)]

                # sort, split, find mean relative errors, and append:
                tup_sorted = sorted(inf_tup, key=lambda tup: tup[1])
                low, high = split_list(tup_sorted)
                guess_influencers, _ = zip(*high)
                guess_others, _ = zip(*low)

                # finde relative errors:
                errors_influencers = RE(d, np.array(guess_influencers))
                errors_noninfluencers = RE(d, np.array(guess_others))

                MRE_i = np.mean(errors_influencers)
                MRE_o = np.mean(errors_noninfluencers)
                mi.append(MRE_i)
                mo.append(MRE_o)

                # find CIs
                CI_influencers = sms.DescrStatsW(errors_influencers).tconfint_mean()
                CI_noninfluencers = sms.DescrStatsW(errors_noninfluencers).tconfint_mean()
                MREs_i_CI.append(CI_influencers)
                MREs_o_CI.append(CI_noninfluencers)

            # after simul repetitions, find mean values and append:
            mean_mi, mean_mo = np.mean(mi), np.mean(mo)
            MREs.append((mean_mi, mean_mo))

            # after simul repetitions, find the 97,5 % percentile and
            # the 2,5 % percentile of mi and mo, which define confidence
            # intervals around the mean (mean_mi and mean_mo):
            # MREs_i_CI.append((np.percentile(mi,2.5), np.percentile(mi,97.5)))
            # MREs_o_CI.append((np.percentile(mo,2.5), np.percentile(mo,97.5)))


    # # unpack values related to high-influencers
    # error_inf, _ = zip(*MREs)
    # inf_lower, inf_upper = zip(*MREs_i_CI)
    # print('high-influencers:', error_inf, inf_lower, inf_upper)
    #
    # # unpack values related to low-influencers
    # _, error_ninf = zip(*MREs)
    # ninf_lower, ninf_upper = zip(*MREs_o_CI)
    # print('low-influencers:', error_ninf, ninf_lower, ninf_upper)

    print(MREs_i_CI)
    # unpack and format confidence intervals
    error_inf, error_ninf = zip(*MREs)
    inf_lower, inf_upper = zip(*MREs_i_CI)
    inf_lower = np.array(error_inf) - np.array(inf_lower)
    inf_upper = np.array(inf_upper) - np.array(error_inf)
    ninf_lower, ninf_upper = zip(*MREs_o_CI)
    ninf_lower = np.array(error_ninf) - np.array(ninf_lower)
    ninf_upper = np.array(ninf_upper) - np.array(error_ninf)
    print('high-influencers:', error_inf, inf_lower, inf_upper)
    print('low-influencers:', error_ninf, ninf_lower, ninf_upper)
    # plot
    x = np.arange(10.0)
    ax.errorbar(x, error_inf, yerr=[inf_lower,inf_upper], capsize=5, fmt='o', c=colors[1], label='high-influencers')
    ax.errorbar(x, error_ninf, yerr=[ninf_lower,ninf_upper], capsize=5, fmt='o',c=colors[0], label='low-influencers')

    x_axis_labels = ['55 dots\n3 views', '55 dots\n9 views',
                     '148 dots\n3 views', '148 dots\n9 views',
                     '403 dots\n3 views', '403 dots\n9 views',
                     '1097 dots\n3 views', '1097 dots\n9 views',
                     '1233 kilo\n3 views', '1233 kilo\n9 views']
    plt.xticks([i for i in range(10)], x_axis_labels, rotation=90)
    ax.yaxis.set_major_formatter(FuncFormatter('{0:.0%}'.format))
    plt.ylabel('Mean relative percentage error')
    legend = ax.legend(title='Simulation results', loc='upper left')
    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    # plt.savefig(os.path.join(PLOTS_DIR, 'figS9.png'), transparent=True, bbox_inches='tight', dpi=300)
    # plt.savefig(os.path.join(PLOTS_DIR, 'figS9.pdf'), transparent=True, bbox_inches='tight', dpi=300)

    plt.show()
    print('average difference btw influencers and non-influencers =',
          np.mean(np.array(error_ninf)-np.array(error_inf)))

# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)
df_all = df_all[df_all.method == 'history']

generate_simulation_data(df_all)
