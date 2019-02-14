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

"""
model simulation + calculation of social influence score:
subjects have an initial "hunch" (taken from the control treatment),
see v previous guesses, take the mean of all numbers and report it as
their estimate.
"""

datafiles = [
            '../data/dots/all_dots_trimmed_anonymous.csv',
            '../data/ox/all_ox_untrimmed_anonymous.csv',
            ]

random.seed(4)
tlength = 400
noise = 0.01
simuls = 1000


def social_influence_scores(guess, seen_guesses):
    influences = np.array([1/(abs(guess-sg) + noise) for sg in seen_guesses])
    normalizer = sum(influences)
    scores = influences/normalizer
    return scores


# mean relative error
def MRE(targets, predictions):
    predictions = np.array(predictions)
    return np.mean(abs((targets - predictions)/targets))


def split_list(sorted_list):
    half = int(len(sorted_list)/2)
    return sorted_list[:half], sorted_list[half:]


def find_social_influence(thread, histories):
    # initialize dict for social influence
    si = dict.fromkeys(list(set(thread)), 0)

    for idx_in_thread, guess in enumerate(thread):
        seen_guesses = histories[idx_in_thread]

        si_score = social_influence_scores(guess, seen_guesses)

        # update the social influence values of the seen guesses:
        for idx, sg in enumerate(seen_guesses):
            si[sg] += si_score[idx]

    # make a new list for the social influence of all guesses
    # Does it? thread[-1] is a guess, so are we automatically setting a guess to have 0 influence?
    si[thread[-1]] = 0

    influence = np.ones(len(thread))

    for idx, guess in enumerate(thread):
        influence[idx] = si[guess]

    return influence


def simulate_thread(views, sampling_pool, thread_length):
    random_draws = random.sample(sampling_pool, thread_length)
    thread = []
    histories = []
    thread.extend(random_draws[:views])  # fill with the first v samples
    histories.extend([random_draws[:views] for i in range(views)])  # fill with the first v histories

    # simulate a thread of length thread_length
    for g in range(views, thread_length):
        history = thread[-views:]  # these are the previous estimates
        own_hunch = random_draws[g]
        final_guess = 1 / (views + 1) * (own_hunch + sum(history))
        thread.append(final_guess)
        histories.append(history)

    return thread, histories


def generate_and_plot(df_all):
    fig, ax = plt.subplots(1,1)
    colors = ['#d94801', '#6a51a3']

    MREs = []
    MREs_i_CI = []
    MREs_o_CI = []
    for position, d in enumerate([55,148,403,1097,1233]):
        # define the control group
        df = df_all[df_all['dots'] == d]
        df_control = df[df.views == 0]
        sampling_pool = list(df_control['guess'].values)

        # print('\n', d, bs.bootstrap(np.array(control), stat_func=bs_stats.mean))

        for v in [3, 9]:
            print("simulating dots: {}, views: {}".format(d, v))
            mi = []
            mo = []

            # simluate with 'simul' runs:
            for i in range(simuls):
                thread, histories = simulate_thread(v, sampling_pool, tlength)
                influence = find_social_influence(thread, histories)

                # make tuples of guesses and their influences
                inf_tup = [(i, influence[pos]) for pos, i in enumerate(thread)]

                # sort, split, find mean relative errors, and append:
                tup_sorted = sorted(inf_tup, key=lambda tup: tup[1])
                low, high = split_list(tup_sorted)
                guess_influencers, _ = zip(*high)
                guess_others, _ = zip(*low)
                MRE_i = MRE(d, np.array(guess_influencers))
                MRE_o = MRE(d, np.array(guess_others))
                mi.append(MRE_i)
                mo.append(MRE_o)


            # after simul repetitions, find mean values and append:
            mean_mi, mean_mo = np.mean(mi), np.mean(mo)
            MREs.append((mean_mi, mean_mo))

            # after simul repetitions, find the 97,5 % percentile and
            # the 2,5 % percentile of mi and mo, which define confidence
            # intervals around the mean (mean_mi and mean_mo):
            MREs_i_CI.append((np.percentile(mi,2.5), np.percentile(mi,97.5)))
            MREs_o_CI.append((np.percentile(mo,2.5), np.percentile(mo,97.5)))

    # unpack values related to high-influencers
    error_inf, _ = zip(*MREs)
    inf_lower, inf_upper = zip(*MREs_i_CI)

    # unpack values related to low-influencers
    _, error_ninf = zip(*MREs)
    ninf_lower, ninf_upper = zip(*MREs_o_CI)

    # plot
    x = np.arange(10.0)
    ax.errorbar(x, error_inf, yerr=[inf_lower,inf_upper], capsize=5, fmt='o', c=colors[1], label='influencers')
    ax.errorbar(x, error_ninf, yerr=[ninf_lower,ninf_upper], capsize=5, fmt='o',c=colors[0], label='non-influencers')

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
    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig('../plots/FigS11.png', transparent=True, dpi=300)
    plt.savefig('../plots/FigS11.pdf', transparent=True, dpi=300)
    plt.show()
    print('average difference btw influencers and non-influencers =',
          np.mean(np.array(error_ninf)-np.array(error_inf)))

# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)
df_all = df_all[df_all['method'] == 'history']

generate_and_plot(df_all)
