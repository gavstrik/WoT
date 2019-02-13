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

file = pd.read_csv('../data/dots/5du4txa7_untrimmed_anonymous.csv') # dots example

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
    player_ids = df['id'].values
    value_nearest_guess = []  # value_nearest_guess[idx] is the value of the guess closest to guess nr. idx
    idx_of_nearest_guess = []  # idx_of_nearest_guess[idx] is the index of the guess closest to guess nr. idx

    for player_id in player_ids:
        own_order = df[df['id'] == player_id]['order'].item() - 1
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
            idx_of_nearest_guess.append((own_order, own_order))
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
        idx_of_nearest_guess.append((own_order, followed_player_order))

    # make a new list for the social influence score of all seen ids
    influence = []
    for player_id in player_ids:
        if player_id in si:
            influence.append(si[player_id])
        else:
            influence.append(0)

    return influence, value_nearest_guess, idx_of_nearest_guess, si_followed


def plotting_tree(df):
    influence, value_nearest_guess, idx_of_nearest_guess, si_followed = calculate_spanning_tree(df)
    guesses = df['guess'].values

    # plot and color initializations
    plt.figure(figsize=(2, 2))
    color_map = plt.get_cmap('magma_r')
    c_norm = colors.Normalize(vmin=0, vmax=max(influence))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=color_map)

    # plot lines between closest guesses, colored according to their influence
    player_ids = df['id'].values
    for player_idx, player_id in enumerate(player_ids):
        plt.plot(
            value_nearest_guess[player_idx],
            idx_of_nearest_guess[player_idx],
            zorder=-1,
            linewidth=3.0,
            color=scalar_map.to_rgba(si_followed[player_idx]))

    # plot the guesses colored with their total social influence score:
    plt.scatter(guesses, [i for i in range(len(guesses))], s=20, c=influence, cmap=color_map)
    plt.xlim(100, 600)
    plt.ylim(12.5, 19.5)

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    output_path = os.path.join(PLOTS_DIR, 'fig4_excerpt.png')

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig(output_path, transparent=True, bbox_inches='tight', dpi=300)
    plt.savefig(output_path, transparent=True, bbox_inches='tight', dpi=300)
    plt.show()

    # Useful for annotating the excerpt in photoshop/illustrator etc.
    # for i in range(13,20):
    #    print(guesses[i], near[i], ords[i], influence[i], si_followed[i])


# main code
df = pd.DataFrame(file)
plotting_tree(df)
