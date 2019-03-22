import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'

"""
This python script plots a spanning tree / hamilton tree from a WoT session
"""

datafiles = [
            '../data/dots.xls',
            '../data/ox.xls',
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
    return idx, array[idx]


def calculate_spanning_tree(df):
    si = {}  # initialize dict for the total social influence score for each player id
    si_followed = []  # array for the social influence of the nearest precursor by position of guess in thread
    player_ids = df.id.values
    value_nearest_guess = []  # value_nearest_guess[idx] is the value of the guess closest to guess nr. idx
    idx_of_nearest_guess = []  # idx_of_nearest_guess[idx] is the index of the guess closest to guess nr. idx

    for guess_idx, player_id in enumerate(player_ids):
        own_guess = df[df['id'] == player_id]['guess'].item()
        seen_ids = df[df['id'] == player_id]['hist'].item()
        # seen_ids is stored in otree as a string, so turn it into a list of player_ids
        seen_ids = pd.eval(seen_ids)
        # Disregard guess in the spanning tree if the data has been trimmed to exclude outliers
        seen_ids = [id for id in seen_ids if id in player_ids]
        seen_guesses = [df[df['id'] == player_id]['guess'].item() for player_id in seen_ids]

        if not seen_guesses:
            # Special case for the first player, who haven't seen any other guesses
            si_followed.append(0)
            value_nearest_guess.append((own_guess, own_guess))
            idx_of_nearest_guess.append((guess_idx, guess_idx))
            continue
        else:
            si_followed.append(max(social_influence_scores(own_guess, seen_guesses)))

        # find the nearest guess
        (nearest_guess_idx, nearest_guess_value) = find_nearest(seen_guesses, own_guess)

        # update the social influence values of the seen ids:
        for idx, player_id in enumerate(seen_ids):
            if player_id not in si:
                si[player_id] = social_influence_scores(own_guess, seen_guesses)[idx]
            else:
                si[player_id] += social_influence_scores(own_guess, seen_guesses)[idx]

        followed_player_order = df[df['id'] == seen_ids[nearest_guess_idx]]['order'].item() - 1

        # append tuples of nearest guesses and orders
        value_nearest_guess.append((own_guess, nearest_guess_value))
        idx_of_nearest_guess.append((guess_idx, followed_player_order))

    # make a new list for the social influence score of all seen ids
    influence = []
    for player_id in player_ids:
        if player_id in si:
            influence.append(si[player_id])
        else:
            influence.append(0)

    return influence, value_nearest_guess, idx_of_nearest_guess, si_followed


def plotting_tree(df_all):
    # plot and color initializations
    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(8.7,12.3))
    ax = ax.flatten()
    color_map = plt.get_cmap('magma_r')
    agg_colors = ["#34495e", '#6a51a3', '#d94801']

    # list the threads
    threads = df_all.session.unique()
    im = [0,0]
    for thread_idx, thread in enumerate(threads):
        df = df_all[df_all.session == thread]
        influence, value_nearest_guess, idx_of_nearest_guess, si_followed = calculate_spanning_tree(df)
        guesses = df.guess.values
        true_number_of_dots = df.d.unique().item()
        views = df.v.unique().item()
        c_norm = colors.Normalize(vmin=0, vmax=max(influence))
        scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=color_map)

        # plot lines between closest guesses, colored according to their influence
        player_ids = df.id.values
        for player_idx, player_id in enumerate(player_ids):
            ax[thread_idx].plot(
                value_nearest_guess[player_idx],
                idx_of_nearest_guess[player_idx],
                zorder=-1, linewidth=2.0,
                color=scalar_map.to_rgba(si_followed[player_idx]))

        # plot the guesses colored with their total social influence score:
        im[thread_idx] = ax[thread_idx].scatter(guesses, [i for i in range(len(guesses))],
                               s=20, c=influence, cmap=color_map)

        # add the moving average and the moving median to the plot
        moving_average = []
        moving_media = []
        for g in range(1, len(guesses)+1):
            guess_list = guesses[:g]
            moving_average.append(np.mean(guess_list))
            moving_media.append(np.median(guess_list))

        ax[thread_idx].plot(moving_average, [i for i in range(len(guesses))], linewidth=1, c=agg_colors[1])
        ax[thread_idx].plot(moving_media, [i for i in range(len(guesses))], linewidth=1, c=agg_colors[2])
        # ax[thread_idx].set_colorbar(shrink=0.5)
        ax[thread_idx].set_xlim(0, 2500)
        #plt.ylim(0, 460)
        ax[thread_idx].axvline(x=true_number_of_dots, linewidth=1, color=agg_colors[0],
                               label='d='+str(true_number_of_dots)+', v='+str(views))
        ax[thread_idx].set_xlabel('estimate')
        ax[thread_idx].legend(loc="upper right", prop={'size': 9})

    fig.colorbar(im[0], ax=ax[0])
    fig.colorbar(im[1], ax=ax[1])
    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig(os.path.join(PLOTS_DIR, 'figS7.png'), transparent=True, bbox_inches='tight', dpi=400)
    plt.savefig(os.path.join(PLOTS_DIR, 'figS7.pdf'), transparent=True, bbox_inches='tight', dpi=400)
    plt.show()


# start here
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_excel(datafile))
    df_all = df_all.append(df, sort=True)

# choose the threads:
df_all = df_all[(df_all.session == 'wv4xujg7') | (df_all.session == 'no1x9itu')]

plotting_tree(df_all)
