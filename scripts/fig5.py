import os
import numpy as np
import pandas as pd
import statsmodels.stats.proportion as smsp
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'

"""
Plotting herding and bandwagoning.
"""

datafiles = [
            '../data/dots.csv',
            '../data/ox.csv',
            ]


def bandwagoning(dfd, p):
    sessions = dfd.session.unique()
    bandwagon_types = {}
    for session in sessions:
        df = dfd[dfd.session == session]
        ids = df['id'].values
        for id in ids:
            own_guess = df[df['id'] == id]['guess'].item()
            seen_ids = df[df['id'] == id]['hist'].item()
            seen_ids = pd.eval(seen_ids)
            seen_ids = [h for h in seen_ids if h in ids]
            if not seen_ids:
                continue
            seen_guesses = [df[df['id'] == g]['guess'].item() for g in seen_ids]

            # bandwagoning: get the unique estimates and their wagon length:
            unique_guesses, counts_unique_guesses = np.unique(seen_guesses, return_counts=True)

            # get the unique wagons:
            unique_wagons = np.unique(counts_unique_guesses)

            # register when a participant jumps on a bandwagon
            length_of_bandwagon_joined = 0
            bonus = df[df['id'] == id]['bonus'].item()

            for idx, guess in enumerate(unique_guesses):
                if guess-(guess*p/100) <= own_guess <= guess+(guess*p/100):
                    length_of_bandwagon_joined = counts_unique_guesses[idx]

            # update the bandwagon dict:
            for length_of_wagon in unique_wagons:
                (times_seen, times_joined, bonus_for_joiners) = bandwagon_types.get(length_of_wagon, (0, 0, 0))
                times_seen += 1

                if length_of_wagon == length_of_bandwagon_joined:
                    times_joined += 1
                    bonus_for_joiners += bonus  # total bonus to bandwagoneers

                bandwagon_types[length_of_wagon] = (times_seen, times_joined, bonus_for_joiners)

    return bandwagon_types


def plot_herders(df_all):
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(11.4,6))
    ax = ax.flatten()
    colors = ['#bcbddc', '#6a51a3', '#d94801']
    views = 9

    # plot cummulative herding as a function of percentage distance
    p_max = 10
    for idx_d, d in enumerate([403, 1097, 1233]):
        df_all = df_all[df_all.v == views]
        dfd = df_all[df_all.d == d]
        N = len(dfd)
        total_bonus = dfd.bonus.sum()
        avg_winrate = total_bonus/N
        win_rate = []
        win_rate_CI = []
        total_herders = []
        total_herders_CI = []
        for p in range(p_max + 1):
            bandwagons = bandwagoning(dfd, p)

            # plot bandwagoning for p = 0:
            if p == 0:
                rate_of_herding = []
                rate_of_herding_CI = []
                for bandwagon_length, (times_seen, times_joined, _) in bandwagons.items():
                    # rate of herding is just the fraction of joiners
                    rate_of_herding.append(times_joined/times_seen)

                    # calculate confidence intervals of binary data, see
                    # https://sigmazone.com/binomial-confidence-intervals/, e.g.
                    # use the "beta"-option giving the Clopper-Pearson exact CI
                    rate_of_herding_CI.append(smsp.proportion_confint(
                                              times_joined, times_seen,
                                              alpha=0.05, method='beta'))

                rate_of_herding_lower, rate_of_herding_upper = zip(*rate_of_herding_CI)
                rate_of_herding_lower = np.array(rate_of_herding) - np.array(rate_of_herding_lower)
                rate_of_herding_upper = np.array(rate_of_herding_upper) - np.array(rate_of_herding)
                x = np.arange(float(len(rate_of_herding)))
                ax[2].errorbar(x, rate_of_herding, yerr=[rate_of_herding_lower,
                                                         rate_of_herding_upper],
                               elinewidth=2, markeredgewidth=2, markersize=8,
                               capsize=5, fmt='x', c=colors[idx_d], label='d='+str(d))

            # get the cumulative rate of herding until p_max:
            number_of_herders = np.sum([times_joined for bandwagon_length, (times_seen, times_joined, _) in bandwagons.items()])
            total_herders.append(number_of_herders/len(dfd))

            # get the number of participant who have seen one or more previous guesses:
            bandwagons_seen = bandwagons[1][0]  # excluding the first estimates of each thread

            # get the confidence intervals, using the normal distribution
            total_herders_CI.append(smsp.proportion_confint(
                                      number_of_herders, bandwagons_seen,
                                      alpha=0.05, method='normal'))

            # get the win rate
            wins = np.sum([bonus_for_joiners for (_, _, bonus_for_joiners) in bandwagons.values()])
            win_rate.append(wins/number_of_herders)
            # get the confidence intervals of the win rate:
            win_rate_CI.append(smsp.proportion_confint(
                                      wins, number_of_herders,
                                      alpha=0.05, method='normal'))

        # plot the cumulative rate of herding including the CIs:
        total_herders_lower, total_herders_upper = zip(*total_herders_CI)
        total_herders_lower = np.array(total_herders) - np.array(total_herders_lower)
        total_herders_upper = np.array(total_herders_upper) - np.array(total_herders)
        x = np.arange(float(len(total_herders)))
        ax[0].errorbar(x, total_herders, yerr=[total_herders_lower,
                                               total_herders_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8,
                       capsize=5, fmt='o', c=colors[idx_d], label='d='+str(d)+' (N='+str(N)+')')

        # plot the relative win rate on the second subplot:
        total_win_rate_lower, total_win_rate_upper = zip(*win_rate_CI)
        total_win_rate_lower = np.array(win_rate) - np.array(total_win_rate_lower)
        total_win_rate_upper = np.array(total_win_rate_upper) - np.array(win_rate)
        x = np.arange(float(len(win_rate)))
        rel_win_rate = win_rate/avg_winrate
        ax[1].errorbar(x, rel_win_rate, yerr=[total_win_rate_lower,
                                               total_win_rate_upper],
                       elinewidth=2, markeredgewidth=2, markersize=8,
                       capsize=5, fmt='v', c=colors[idx_d], label='d='+str(d)+' (N='+str(N)+')')
        ax[1].plot(x, np.poly1d(np.polyfit(x, rel_win_rate, 1))(x), color=colors[idx_d])

    # plotting paraphernalia
    ax[0].legend(loc='lower right', prop={'size': 11})
    ax[0].set_xlabel('p (%)', fontsize='x-large')
    ax[0].set_ylabel('Fraction of herders', fontsize='x-large')
    ax[0].set_xticks([0, 2, 4, 6, 8, 10])
    ax[0].set_xticklabels(['0', '2', '4', '6', '8', '10'], fontsize='large')
    ax[0].set_yticks([i/5 for i in range(6)])
    ax[0].set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'], fontsize='large')
    ax[0].set_ylim(-0.03,1.03)

    ax[1].legend(loc='lower right', prop={'size': 11})
    ax[1].set_xlabel('p (%)', fontsize='x-large')
    ax[1].set_ylabel('Relative win rate', fontsize='xx-large')
    ax[1].set_xticks([0, 2, 4, 6, 8, 10])
    ax[1].set_xticklabels(['0', '2', '4', '6', '8', '10'], fontsize='large')
    ax[1].set_yticks([.6, .8, 1., 1.2, 1.4, 1.6])
    ax[1].set_yticklabels(['0.6', '0.8', '1.0', '1.2', '1.4', '1.6'], fontsize='large')
    ax[1].set_ylim(0.57,1.63)

    ax[2].legend(loc='upper left', bbox_to_anchor=(0.0,0.903), prop={'size': 11})
    ax[2].set_xlabel('# identical estimates (p = 0)', fontsize='x-large')
    ax[2].set_ylabel('Rate of herding', fontsize='xx-large')
    ax[2].set_xticks([i for i in range(len(bandwagons))])
    ax[2].set_xticklabels(['1', '2', '3', '4', '5', '6', '7', '8', '9'], fontsize='large')
    ax[2].set_yticks([i/5 for i in range(6)])
    ax[2].set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'], fontsize='large')
    ax[2].set_ylim(-0.03,1.03)

    txtA = fig.text(0.08, .904, 'A', fontsize='xx-large', fontweight='bold')
    txtB = fig.text(0.405, .904, 'B', fontsize='xx-large', fontweight='bold')
    txtB = fig.text(0.73, .904, 'C', fontsize='xx-large', fontweight='bold')
    plt.grid(False)
    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig(os.path.join(PLOTS_DIR, 'fig5.png'), transparent=True, bbox_extra_artists=(txtA,txtB),
                bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(PLOTS_DIR, 'fig5.pdf'), transparent=True, bbox_extra_artists=(txtA,txtB),
                bbox_inches='tight', dpi=300)
    plt.show()


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.read_csv(datafile)
    df_all = df_all.append(df, sort=True)
df_all = df_all[df_all.method == 'history']

plot_herders(df_all)
