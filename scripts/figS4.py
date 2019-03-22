import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams["font.family"] = "sans-serif"
PLOTS_DIR = '../plots'


"""
Plotting the probability of being better than the mean of estimates seen.
"""

datafiles = [
            '../data/dots.xls',
            '../data/ox.xls',
            ]


# relative error
def RE(targets, predictions):
    return abs(targets - predictions)/targets


# mean relative error
def MRE(targets, predictions):
    predictions = np.array(predictions)
    return np.mean(abs(targets - predictions)/targets)


# median relative error
def MdRE(targets, predictions):
    predictions = np.array(predictions)
    return np.median(abs(targets - predictions)/targets)


# outliers with an error rate above 10 are removed
def remove_outliers(df, true_number_of_dots):
    max_error_rate = 10
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    df = df[(df.guess >= lower) & (df.guess <= upper)]
    return df


def better_than_sample_median(df_all):
    betterBars = []
    equalBars = []
    worseBars = []
    how_much_better_line_mean = []
    how_much_better_line_median = []
    for d in [55, 148, 403, 1097, 1233]:
        better_than_median = []
        better_than_mean = []
        how_much_better_than_median = []
        how_much_better_than_mean = []
        df = df_all[df_all.d == d]
        df = remove_outliers(df, d)
        participants = df.code.values
        within = 0


        for participant in df.code.values:
            own_guess = df[df.code == participant]['guess'].item()
            if own_guess <= 0:
                continue
            seen_ids = pd.eval(df[df.code == participant]['hist'].item())
            if not seen_ids:
                continue
            session = df[df.code == participant]['session'].item()
            seen_participants = []
            for g in seen_ids:
                if not df[(df.session == session) & (df.id == g)]['code'].values:
                    continue
                seen_participants.append(df[(df.session == session) & (df.id == g)]['code'].item())

            # remove id's from hist which have been trimmed away:
            seen_participants = [h for h in seen_participants if h in participants]

            if not seen_participants:
                continue
            seen_guesses = [df[df.code == g]['guess'].item() for g in seen_participants]
            error_guess = RE(d, own_guess)

            # check how many estimates are within their sample limits:
            if own_guess >= np.min(seen_guesses) and own_guess <= np.max(seen_guesses):
                within += 1

            # check if error of estimate is smaller than the median sample error:
            median_sample_error = MdRE(d, seen_guesses)
            if error_guess < median_sample_error:
                better = 1
            elif error_guess == median_sample_error:
                better = 0
            else:
                better = -1
            better_than_median.append(better)
            how_much_better_than_median.append(median_sample_error - error_guess)

            # check if error of estimate is smaller than the mean sample error:
            mean_sample_error = MRE(d, seen_guesses)
            if error_guess < mean_sample_error:
                better = 1
            elif error_guess == mean_sample_error:
                better = 0
            else:
                better = -1
            better_than_mean.append(better)
            how_much_better_than_mean.append(mean_sample_error - error_guess)

        # print(d, within, len(participants), within/len(participants))

        # only use means for plotting:
        betterBars.append(sum([1 for k in better_than_mean if k == 1]))
        equalBars.append(sum([1 for k in better_than_mean if k == 0]))
        worseBars.append(sum([1 for k in better_than_mean if k == -1]))

        # append the mean and median improvement of subsequent guesses
        how_much_better_line_mean.append(np.mean(how_much_better_than_mean))
        how_much_better_line_median.append(np.mean(how_much_better_than_median))

    # plotting stuff
    r = [0,1,2,3,4]  # "x-axis"
    raw_data = {'betterBars': betterBars, 'equalBars': equalBars, 'worseBars': worseBars}
    df_plot = pd.DataFrame(raw_data)
    colors = ["#d94801", "#95a5a6", "#807dba", '#5ab4ac', "#d8b365"]

    # From raw value to percentage
    totals = [i+j+k for i, j, k in zip(df_plot['betterBars'], df_plot['equalBars'],
                                       df_plot['worseBars'])]
    betterBars = [i/j for i, j in zip(df_plot['betterBars'], totals)]
    equalBars = [i/j for i, j in zip(df_plot['equalBars'], totals)]
    worseBars = [i/j for i, j in zip(df_plot['worseBars'], totals)]

    # plot
    barWidth = 0.85
    names = ('55 dots', '148 dots', '403 dots', '1097 dots', '1233 kilo')
    fig, ax = plt.subplots(1,1)

    # make barplots
    b1 = ax.bar(r, worseBars, bottom=[i+j for i,j in zip(betterBars, equalBars)],
            color=colors[0], edgecolor='white', width=barWidth, label='smaller')
    b2= ax.bar(r, equalBars, bottom=betterBars, color=colors[1], edgecolor='white',
            width=barWidth, label='same')
    b3 = ax.bar(r, betterBars, color=colors[2], edgecolor='white', width=barWidth,
                label='higher')
    plt.axhline(y=.5, xmax=0.95, linewidth=1, color='black')
    plt.xticks(r, names)
    ax.set_yticks([-0.05, 0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1])
    ax.set_yticklabels(['', '0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%',
                        '80%', '90%', '100%'])
    plt.ylabel('Estimation error vs. mean sample error (RE)')

    # plot the improvements
    ax2 = ax.twinx()
    #b4 = ax2.plot(how_much_better_line_median, linestyle='-', marker='o',
    #              label='median', c=colors[3])
    b5 = ax2.plot(how_much_better_line_mean, linestyle='-', marker='o',
                  label='mean', c=colors[4])
    ax2.axhline(y=0, xmin=0.05, linestyle='--', linewidth=1, color='black')
    ax2.set_yticks([-.01, 0, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1])
    ax2.set_yticklabels(['', '0%', '1%', '2%', '3%', '', '', '', '', '', '', ''])
    ax2.set_ylabel('Mean improvement\nof next estimate (RE)')
    ax2.yaxis.set_label_coords(1.1,.23)
    ax.grid(False), ax2.grid(False), plt.grid(False)

    # control position of labels by sorting both labels and handles by labels
    handles, labels = ax.get_legend_handles_labels()
    handles_lines, labels_lines = ax2.get_legend_handles_labels()
    handles.extend(handles_lines)
    labels.extend(labels_lines)
    plt.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.tight_layout()

    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig(os.path.join(PLOTS_DIR, 'figS4.png'), transparent=True, dpi=300)
    plt.savefig(os.path.join(PLOTS_DIR, 'figS4.pdf'), transparent=True, dpi=300)
    plt.show()


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.read_excel(datafile)
    df_all = df_all.append(df, sort=True)
df_all = df_all[df_all.method == 'history']

df_all = df_all[df_all.v != 0]
better_than_sample_median(df_all)
