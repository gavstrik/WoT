import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
# plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"

"""
Plotting herding and bandwagoning.
"""

datafiles = [
            '../data/dots/all_dots_untrimmed_anonymous.csv',
            '../data/ox/all_ox_untrimmed_anonymous.csv',
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
            # print(unique_wagons)

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
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7,5))
    ax = ax.flatten()
    colors = ['#bcbddc', '#6a51a3', '#d94801']
    views = 9

    # plot cummulative herding as a function of percentage distance
    p_max = 10
    for pos, d in enumerate([403, 1097, 1233]):
        df_all = df_all[df_all.views == views]
        dfd = df_all[df_all.dots == d]
        total_bonus = dfd['bonus'].sum()
        avg_winrate = total_bonus/len(dfd)

        total_herders = []
        win_rate = []
        for p in range(p_max + 1):
            bandwagons = bandwagoning(dfd, p)

            # plot bandwagoning for p = 0:
            if p == 0:
                band = []
                for bandwagon_length, (times_seen, times_joined, _) in bandwagons.items():
                    band.append(times_joined/times_seen)

                ax[1].plot(band, color=colors[pos], label='d = '+str(d))

            # get the cumulative rate of herding until p_max:
            herders = np.sum([times_joined for bandwagon_length, (times_seen, times_joined, _) in bandwagons.items()])
            total_herders.append(herders/len(dfd))

            # get the win rate for d = 1233
            if d == 1233:
                wins = np.sum([bonus_for_joiners for (_, _, bonus_for_joiners) in bandwagons.values()])
                win_rate.append(wins/herders)

        # plot the cumulative rate of herding
        ax[0].plot(total_herders, color=colors[pos], label='d = '+str(d))

    # plot the relative win rate for d = 1233 on a secondary y-axis:
    rel_win_rate = win_rate/avg_winrate
    # print('relative win rate of herders:\n', rel_win_rate)
    ax0 = ax[0].twinx()
    ax0.plot(rel_win_rate, color=colors[pos], linestyle='--', label='rel. win rate')
    ax0.set_yticks([1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45])

    # plotting paraphernalia
    handles, labels = ax[0].get_legend_handles_labels()
    handles_winrate, labels_winrate = ax0.get_legend_handles_labels()
    handles.extend(handles_winrate)
    labels.extend(labels_winrate)
    ax[0].legend(handles, labels, loc='lower right', ncol=1)
    ax[0].set_xlabel('p', fontsize='large')
    ax[0].set_ylabel('Fraction of herders', fontsize='large')

    ax[1].set_xticks([i for i in range(len(bandwagons))])
    ax[1].set_xticklabels(['1', '', '3', '', '5', '', '7', '', '9'])
    ax[1].legend(loc='lower right') #, prop={'size': 8})
    ax[1].set_xlabel('Identical estimates (p = 0)', fontsize='large')
    ax[1].set_ylabel('Rate of herding', fontsize='large')
    txtA = fig.text(0.12, .904, 'A', fontsize='xx-large', fontweight='bold')
    txtB = fig.text(0.64, .904, 'B', fontsize='xx-large', fontweight='bold')
    plt.grid(False)
    plt.tight_layout()
    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig('../plots/Fig6.png', transparent=True, bbox_extra_artists=(txtA,txtB),
                bbox_inches='tight', dpi=300)
    plt.savefig('../plots/Fig6.pdf', transparent=True, bbox_extra_artists=(txtA,txtB),
                bbox_inches='tight', dpi=300)
    plt.show()


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.read_csv(datafile)
    df_all = df_all.append(df)
df_all = df_all[df_all.method == 'history']

plot_herders(df_all)
