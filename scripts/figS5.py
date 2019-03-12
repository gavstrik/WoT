import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'


"""
Histogram of the differences between individual estimates and sample
means/medians showing that most estimates are close to the mean and
the median of the sample seen.
"""

datafiles = [
            '../data/dots/all_dots_untrimmed_anonymous.csv',
            '../data/ox/all_ox_untrimmed_anonymous.csv',
            ]


def error_rate(prediction, reference):
    return abs(prediction-reference)/reference


# making histogram of variances independent of image
def histograms(df):
    median_diff = []
    mean_diff = []
    colors = ['#6a51a3', '#d94801']

    codes = df.code.values  # these are all the individuals having made a guess
    within = 0
    for code in codes:
        own_guess = df[df.code == code]['guess'].item()
        seen_ids = pd.eval(df[df.code == code]['hist'].item())
        session = df[df.code == code]['session'].item()

        if not seen_ids:
            continue

        seen_codes = [df[(df.session == session) & (df.id == g)]['code'].item() for g in seen_ids]
        seen_guesses = [df[(df.code == g) & (df.session == session)]['guess'].item() for g in seen_codes]
        mean_diff.append(error_rate(own_guess, np.mean(seen_guesses)))
        median_diff.append(error_rate(own_guess, np.median(seen_guesses)))

        if np.min(seen_guesses) <= own_guess <= np.max(seen_guesses):
            within += 1

    # plot
    bins = np.linspace(0, 2, 30)
    weights_mean = np.ones_like(mean_diff)/float(len(mean_diff))
    weights_median = np.ones_like(median_diff)/float(len(median_diff))
    plt.hist([mean_diff, median_diff], bins=bins, stacked=False,
             weights=[weights_mean, weights_median],
             color=colors, label=['sample mean', 'sample median'])

    # plotting paraphernalia
    plt.legend()
    plt.xlabel('error rate btw. estimate and sample mean/median')
    plt.ylabel('fraction of participants')
    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig(os.path.join(PLOTS_DIR, 'figS5.png'), transparent=True, dpi=300)
    plt.savefig(os.path.join(PLOTS_DIR, 'figS5.pdf'), transparent=True, dpi=300)

    plt.show()


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_csv(datafile))
    df_all = df_all.append(df)
df_all = df_all[df_all.method == 'history']
df_all = df_all[df_all.views != 0]
df_all = df_all[df_all.views != 27]

histograms(df_all)
