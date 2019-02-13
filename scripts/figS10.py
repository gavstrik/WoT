import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import matplotlib
# plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"


"""
here we simulate the model: subjects have a 'hunch'
(taken from the control treatment), see v previous guesses, take the mean
of all numbers and report it as their estimate.
"""

datafiles = [
            '../data/dots/all_dots_trimmed_anonymous.csv',
            '../data/ox/all_ox_untrimmed_anonymous.csv',
            ]

tlength = 400
simuls = 1000
random.seed(4)

def means_and_medians(df_t):
    sessions = df_t['session'].unique()
    true_medians = []
    true_means = []
    true_stds = []
    for v in [0,1,3,9]:
        true_medians.append(df_t[df_t.views == v]['guess'].median())
        true_means.append(df_t[df_t.views == v]['guess'].mean())
        true_stds.append(df_t[df_t.views == v]['guess'].std())
    return true_medians, true_means, true_stds


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)
df_all = df_all[df_all['method'] == 'history']


fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(11,6))
axes[-1, -1].axis('off')  # do not show the last subplot
axes = axes.flatten()
colors = ['#6a51a3', '#d94801']

for position, d in enumerate([55, 148, 403, 1097, 1233]):
    # define the control group
    df = df_all[df_all['dots'] == d]
    true_medians, true_means, _ = means_and_medians(df)
    df_control = df[df.views == 0]
    control = df_control['guess'].values

    # simluate with many runs
    medians = []
    medians.append(np.median(control))
    means = []
    means.append(np.mean(control))

    for v in [1, 3, 9]:
        m = []
        md = []

        for i in range(simuls):
            random_draws = random.sample(list(control), tlength)
            thread = []
            thread.extend(random_draws[:v])  # fill thread with the first v samples

            # simulate a thread of length tlength
            for g in range(v, tlength):
                sum_of_seen_history = sum(thread[-v:])  # these are the previous estimates
                own_hunch = random_draws[g]
                final_guess = 1 / (v + 1) * (own_hunch + sum_of_seen_history) # add your own guess
                thread.append(final_guess)

            m.append(np.mean(thread))  # store the new number in the m-thread
            md.append(np.median(thread))  # store the new number in the md-thread

        medians.append(np.mean(md))
        means.append(np.mean(m))

    axes[position].plot(medians, marker='o', linestyle='--', c=colors[1], label='simulated medians')
    axes[position].plot(true_medians, marker='o', linestyle='-', c=colors[1], label='true medians')
    axes[position].plot(means, marker='o', linestyle='--', c=colors[0], label='simulated means')
    axes[position].plot(true_means, marker='o', linestyle='-', c=colors[0], label='true means')

# plotting paraphernalia
handles, labels = axes[0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.72,0.27), title='legend', fontsize='large', ncol=1)
plt.tick_params(axis="x", length=0)  # remove the small ticks
x_axis_labels = ['0', '1', '3', '9']  # rename xticks
plt.xticks([i for i in range(4)], x_axis_labels)

# remove small ticks
for i in range(5):
    axes[i].tick_params(axis="x", length=0)

axes[0].set_title('d = 55', fontsize='large', fontweight='bold')
axes[1].set_title('d = 148', fontsize='large', fontweight='bold')
axes[2].set_title('d = 403', fontsize='large', fontweight='bold')
axes[3].set_title('d = 1097', fontsize='large', fontweight='bold')
axes[4].set_title('d = 1233', fontsize='large', fontweight='bold')
txtv = fig.text(0.54, .0, 'v', fontsize='x-large', ha='center')
plt.tight_layout()
# Remember: save as pdf and transparent=True for Adobe Illustrator
fig.savefig('../plots/FigS10.png', transparent=True, bbox_extra_artists=(txtv,),
            bbox_inches='tight', dpi=300)
plt.show()
