import numpy as np
import pandas as pd
from tabulate import tabulate
from scipy.stats import kurtosis, entropy
from scipy.stats import skew as skewness
from collections import Counter

"""
Pretty printing the most important stats.
Usage: copy output into latex file.
"""


datafiles = [
            '../CleanedData/dots.xlsx',
            ]


def tabulated_stats(df1):
    table = pd.DataFrame()
    dots = []
    method = []
    views = []
    N = []
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
    sessions = df1['session'].unique()
    N_hist = N_max = 0
    turkerdict = {k:0 for k in df1.hashed_turker.unique()}

    for session in sessions:
        df = df1[df1.session == session]
        dots.append(df.d.unique().item())
        guesses = df.guess.values
        views.append(df.v.unique().item())
        method.append(df.method.unique().item())
        N.append(len(df.guess.values))
        median.append(np.around(np.median(guesses),2))
        mean.append(np.around(np.mean(guesses),2))
        sd.append(np.around(np.std(guesses),2))
        cv.append(np.around(np.std(guesses)/np.mean(guesses),2))
        skew.append(np.around(skewness(guesses),2))
        kurt.append(np.around(kurtosis(guesses),2))
        bonus.append(np.around(100*df['bonus'].sum()/len(guesses),2))

        if df.method.unique().item() == 'history':
            N_hist += len(df.guess.values)
        else:
            N_max += len(df.guess.values)

        for turker in df.hashed_turker.unique():
            turkerdict[turker] += 1

    c={}
    for v in turkerdict.values():
        if not v in c:
            c[v]=1
        else:
            c[v]+=1
    print('number of turkers seeing one, two, three or four images:', dict(c.items()))

    table['method'] = method
    table['d'] = dots
    table['v'] = views
    table['thread'] = sessions
    table['N'] = N
    table['median'] = median
    table['mean'] = mean
    table['SD'] = sd
    table['CV'] = cv
    table['skew'] = skew
    table['kurt'] = kurt
    table['bonus (\%)'] = bonus


    table = table.sort_values(['method', 'd', 'v'])
    return table, np.sum(N), N_hist, N_max


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_excel(datafile))
    df_all = df_all.append(df)
table, datapoints, N_hist, N_max = tabulated_stats(df_all)

print('total number of data points:', datapoints)
print('total number of data points in history threads:', N_hist)
print('total number of data points in manipulated threads:', N_max)
print('unique turkers in total:', len(df_all.hashed_turker.unique()))
print(tabulate(table, tablefmt="latex_raw", headers="keys", showindex=False))
