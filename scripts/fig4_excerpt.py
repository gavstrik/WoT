import numpy as np
import pandas as pd
import statistics
from numpy import percentile
from fractions import Fraction
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.pyplot import figure
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import style
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


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
    return (idx, array[idx])


def calculate_spanning_tree(df):
    si = {}  # initialize dict for the total social influence score
    si_followed = []  # array for the social influence of the nearest precusor
    views = df['views'].unique().item()
    ids = df['id'].values
    near = []  # list of tuples containing nearest estimates
    ords = []  # list of tuples containing nearest orders
    for id in ids:
        own_order = df[df['id'] == id]['order'].item() - 1
        own_guess = df[df['id'] == id]['guess'].item()
        seen_ids = df[df['id'] == id]['hist'].item()
        seen_ids = pd.eval(seen_ids)
        seen_ids = [h for h in seen_ids if h in ids]
        seen_guesses = [df[df['id'] == g]['guess'].item() for g in seen_ids]
        if not seen_guesses:
            si_followed.append(0)
            near.append((own_guess, own_guess))
            ords.append((own_order, own_order))
            continue
        else:
            si_followed.append(max(social_influence_scores(own_guess, seen_guesses)))

        # find the nearest guess
        nearest = find_nearest(seen_guesses, own_guess)

        # update the social influence values of the seen ids:
        for pos, g in enumerate(seen_ids):
            if not g in si:
                si[g] = social_influence_scores(own_guess, seen_guesses)[pos]
            else:
                si[g] += social_influence_scores(own_guess, seen_guesses)[pos]

        followed_player_order = df[df['id'] == seen_ids[nearest[0]]]['order'].item() -1

        # append tuples of nearest guesses and orders
        near.append((own_guess, nearest[1]))
        ords.append((own_order, followed_player_order))

    # make a new list for the social influence score of all seen ids
    infl = []
    for id in ids:
        if id in si:
            infl.append(si[id])
        else:
            infl.append(0)

    return infl, near, ords, si_followed


def plotting_tree(df, views, method):
    influence, near, ords, si_followed = calculate_spanning_tree(df)
    print('LEN:', len(df), len(influence), len(near), len(ords), len(si_followed))
    session = df['session'].unique().item()
    guesses = df['guess'].values
    true_number_of_dots = df['dots'].unique().item()

    # plot and color initializations
    plt.figure(figsize=(2, 2))
    colmap = 'magma_r'
    mycmap = cm = plt.get_cmap(colmap)
    cNorm = colors.Normalize(vmin=0, vmax=max(influence))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=mycmap)
    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]

    # plot lines btw nearest guesses, colored according to their influence
    ids = df['id'].values
    print('ids:', len(ids), len(guesses))
    for pos, id in enumerate(ids):
        colorVal = scalarMap.to_rgba(si_followed[pos])
        plt.plot(near[pos], ords[pos], zorder=-1, linewidth=3.0, color=colorVal)

    # plot the guesses colored with their total social influence score:
    plt.scatter(guesses, [i for i in range(len(guesses))], s=30, c=influence, cmap=mycmap)

    # add the moving average and the moving median to the plot
    ma = []
    me = []
    for g in range(1, len(guesses)+1):
        guess_list = guesses[:g]
        ma.append(np.mean(guess_list))
        me.append(np.median(guess_list))
    # plt.plot(ma, [i for i in range(len(guesses))], linewidth=.5, c=flatui[5])
    # plt.plot(me, [i for i in range(len(guesses))], linewidth=.5, c=flatui[3])

    plt.xlim(100, 600)
    plt.ylim(12.5, 19.5)
    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig('../plots/Fig4_excerpt.png', transparent=True, dpi=300)
    plt.show()
    for i in range(13,20):
        print(guesses[i], near[i], ords[i], influence[i], si_followed[i])


# main code
df = pd.DataFrame(file)
plotting_tree(df, df['views'].unique().item(), df['method'].unique().item())
