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
    bandwagons = {}
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
            unique_elements, counts_elements = np.unique(seen_guesses, return_counts=True)

            # get the unique wagons:
            unique_wagons = np.unique(counts_elements)
            # print(unique_wagons)

            # register when a participant jumps on a bandwagon
            bandwagoneer = 0
            bandwagon_size = 0
            bonus = df[df['id'] == id]['bonus'].item()
            for pos, i in enumerate(unique_elements):
                if own_guess >= i-(i*p/100) and own_guess <= i+(i*p/100):
                    bandwagoneer = 1
                    bandwagon_size = counts_elements[pos]

            # update the bandwagon dict:
            for key in unique_wagons:
                (a, b, c) = bandwagons.get(key, (0,0,0))
                a += 1
                if key == bandwagon_size:
                    b += bandwagoneer
                    c += bonus  # total bonus to bandwagoneers
                bandwagons[key] = (a, b, c)
    return bandwagons


def plot_herders(df_all):
    fig,ax = plt.subplots(nrows=1, ncols=2, figsize=(7,5))
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
                bonus = []
                for k, v in bandwagons.items():
                    band.append(v[1]/v[0])
                    # if v[1] != 0:
                    #     bonus.append(v[2]/v[1])
                    # else:
                    #     bonus.append(0)
                ax[1].plot(band, color=colors[pos], label='d = '+str(d))
                # ax1 = ax[1].twinx()
                # bonus = bonus/avg_winrate
                # print(bonus)
                # ax1.plot(bonus, color=colors[pos], linestyle='--', label='d = '+str(d))

            # get the cumulative rate of herding until p_max:
            herders = np.sum([v[1] for k, v in bandwagons.items()])
            total_herders.append(herders/len(dfd))

            # get the win rate for d = 1233
            if d == 1233:
                wins = np.sum([v[2] for k, v in bandwagons.items()])
                win_rate.append(wins/herders)

        # plot the cumulative rate of herding
        ax[0].plot(total_herders, color=colors[pos], label='d = '+str(d))

    # plot the relative win rate for d = 1233 on a secondary y-axis:
    rel_win_rate = win_rate/avg_winrate
    # print('relative win rate of herders:\n', rel_win_rate)
    ax0 = ax[0].twinx()
    ax0.plot(rel_win_rate, color=colors[pos], linestyle='--', label='win rate')
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
