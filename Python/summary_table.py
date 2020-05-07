import numpy as np
import pandas as pd
from tabulate import tabulate
from scipy.stats import kurtosis
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
    skew = []
    kurt = []
    bonus = []
    sessions = df1['session'].unique()
    turkerdict = {k:0 for k in df1.hashed_turker.unique()}

    for session in sessions:
        df = df1[df1.session == session]
        d = df.d.unique().item()
        v = df.v.unique().item()
        dots.append(d)
        guesses = df.guess.values
        views.append(v)
        method.append(df.method.unique().item())
        N.append(len(df.guess.values))
        median.append(np.around(np.median(guesses),2))
        mean.append(np.around(np.mean(guesses),2))
        sd.append(np.around(np.std(guesses),2))
        cv.append(np.around(np.std(guesses)/np.mean(guesses),2))
        skew.append(np.around(skewness(guesses),2))
        kurt.append(np.around(kurtosis(guesses),2))
        tmp = 0
        for g in guesses:
            if (d - d/10) <= g <= (d + d/10):
                tmp += 1
        bonus.append(np.around(100*tmp/len(guesses),2))

        for turker in df.hashed_turker.unique():
            turkerdict[turker] += 1

    c = {}
    for tmp1 in turkerdict.values():
        if not tmp1 in c:
            c[tmp1]=1
        else:
            c[tmp1]+=1
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
    return table


# main code
df_all = pd.DataFrame()
for datafile in datafiles:
    df = pd.DataFrame(pd.read_excel(datafile))
    df_all = df_all.append(df)

print('total number of data points/turkers:',
                len(df_all), len(df_all.hashed_turker.unique()))
print('total number of data points/turkers in history threads:',
                len(df_all[(df.method == 'history') & (df.v != 0)]),
                len(df_all[(df.method == 'history') & (df.v != 0)]['hashed_turker'].unique()))
print('total number of data points/turkers in manipulated threads:',
                len(df_all[df.method == 'max']),
                len(df_all[df.method == 'max']['hashed_turker'].unique()))
print('total number of data points/turkers in control threads:',
                len(df_all[df.v == 0]),
                len(df_all[df.v == 0]['hashed_turker'].unique()))


table = tabulated_stats(df_all)
print(tabulate(table, tablefmt="latex_raw", headers="keys", showindex=False))
