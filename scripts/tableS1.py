import numpy as np
import pandas as pd
from tabulate import tabulate
from scipy.stats import kurtosis, entropy
from scipy.stats import skew as skewness

"""
Pretty printing the most important stats.
Usage: copy output into markdown file.
"""


datafiles = [
            '../data/dots.xls',
            '../data/ox.xls',
            ]


def mae(predictions, truth):
    return np.mean(abs(truth - predictions))/truth


def collective_error_rate(prediction, truth):
    return abs(truth - prediction)/truth


# outliers with an error rate above 10 are removed
def remove_outliers(data, true_number_of_dots):
    max_error_rate = 10
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    newdata = [i for i in data if i >= lower and i <= upper]
    return np.array(newdata)


def tabulated_stats(df1):
    table = pd.DataFrame()
    dots = []
    method = []
    views = []
    length = []
    median = []
    mean = []
    sd = []
    cv = []
    col_error_mean = []
    col_error_median = []
    maes = []
    skew = []
    kurt = []
    bonus = []
    ent = []
    outliers = []
    sessions = df1['session'].unique()
    for session in sessions:
        df = df1[df1.session == session]
        dots.append(df.d.unique().item())
        guesses = remove_outliers(df.guess.values, df.d.unique().item())
        outliers.append(len(df.guess.values)-len(guesses))
        # method.append(df.method.unique().item())
        views.append(df.v.unique().item())
        method.append(df.method.unique().item())
        length.append(len(df.guess.values))
        median.append(np.around(np.median(guesses),2))
        mean.append(np.around(np.mean(guesses),2))
        sd.append(np.around(np.std(guesses),2))
        cv.append(np.around(np.std(guesses)/np.mean(guesses),2))
        col_error_mean.append(np.around(collective_error_rate(np.mean(guesses), df.d.unique().item()),2))
        col_error_median.append(np.around(collective_error_rate(np.median(guesses), df.d.unique().item()),2))
        maes.append(np.around(mae(guesses, df.d.unique().item()),2))
        skew.append(np.around(skewness(guesses),2))
        kurt.append(np.around(kurtosis(guesses),2))
        bonus.append(np.around(100*df['bonus'].sum()/len(guesses),2))
        # ent.append(entropy(guesses, qk=[df.d.unique().item() for i in range(len(guesses))], base=None))

    table['d'] = dots
    table['v'] = views
    table['method'] = method
    table['thread'] = sessions
    table['N'] = length
    table['out'] = outliers
    table['median'] = median
    table['mean'] = mean
    table['err-med'] = col_error_median
    table['err-mean'] = col_error_mean
    table['SD'] = sd
    table['CV'] = cv
    table['MAE'] = maes
    table['skew'] = skew
    table['kurt'] = kurt
    table['bonus (\%)'] = bonus
    # table['KL'] = ent

    table = table.sort_values(['d', 'v', 'thread'])
    return table, np.sum(length), np.sum(outliers)


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_excel(datafile))
    df_all = df_all.append(df)

table, datapoints, outliers = tabulated_stats(df_all)
print(tabulate(table, tablefmt="latex_raw", headers="keys", showindex=False))
print('total number of data points:', datapoints)
print('total number of outliers:', outliers)
