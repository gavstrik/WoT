import numpy as np
from numpy import percentile
import pandas as pd
import hashlib
from collections import Counter

datafiles = [
    '../RawData/all_apps_wide_2018-08-14.csv',
    '../RawData/all_apps_wide_2018-09-13.csv',
    '../RawData/all_apps_wide_2018-10-02.csv',
    '../RawData/all_apps_wide_2019-06-20.csv',
]


# anonymizing turker id with salted hash function:
def hash_simple(df, column):
    return df[column].apply(hash)


def make_dataframe(file):
    df = pd.DataFrame(file)

    # take only the relevant columns in data sheet
    df = df[[
            'participant._current_app_name',
            'session.config.selection_method',
            'session.config.true_number_of_dots',
            'session.config.views',
            'session.code',
            'participant.mturk_worker_id',
            'dots.1.player.decision_order',
            'dots.1.player.decision_history',
            'dots.1.player.guess',
            'dots.1.player.payoff',
            ]]

    # rename columns
    df.columns = [
                  'task',
                  'method',
                  'd',
                  'v',
                  'session',
                  'mturker',
                  'decision_order',
                  'history',
                  'guess',
                  'bonus'
    ]
    # drop all rows with nan's
    df = df.dropna()
    return df

def remove_duplicates(df_all):
    images = df_all.d.unique()
    newdf = pd.DataFrame()

    for image in images:
        df_i = df_all[df_all['d'] == image]
        duplicates = df_i[df_i['mturker'].duplicated(keep=False)]['mturker'].unique()
        df_i = df_i.drop_duplicates(['mturker'], keep='first')
        newdf = newdf.append(df_i)

    print(len(df_all.guess.values) - len(newdf.guess.values), 'duplicates removed...')
    return newdf


def make_cleaned_data(df):
    df_clean = pd.DataFrame()
    sessions = df.session.unique()
    mturker_count = {}
    for session in sessions:
        df_s = df[df['session'] == session]
        mturkers = df_s.mturker.unique()
        for m in mturkers:
            if not m in mturker_count:
                mturker_count[m] = 1
            else:
                mturker_count[m] += 1

        # # sort with "order"
        # df_s = df_s.sort_values('decision_order')

        # save
        df_clean = df_clean.append(df_s)

    print('\ndistribution of mturkers over seen images:')
    print(Counter(mturker_count.values()))

    # make anonymous and drop turker id column
    df_clean['hashed_turker'] = hash_simple(df_clean, 'mturker')
    df_clean = df_clean.drop(['mturker'], 1)

    # sort, reset index, and save
    df_clean = df_clean.sort_values(by=['method', 'd', 'v', 'session', 'decision_order'])
    df_clean = df_clean.reset_index(drop=True)
    df_clean.to_excel('../CleanedData/dots.xlsx')
    print('dots.xlsx saved containing', len(df_clean.guess.values), 'estimates \n \
    from', len(mturker_count.values()), 'turkers.')


# start here
df_all = pd.DataFrame()

for datafile in datafiles:
    df = make_dataframe(pd.read_csv(datafile))
    df_all = df_all.append(df)

# take only the article-relevant runs:
df_all = df_all[(df_all.method == 'history') | (df_all.method == 'max')]
df_all = df_all[df_all.v != 27]
print(len(df_all.guess.values), 'estimates found and...\n',
      len(df_all.mturker.unique()), 'unique respondents found...')

# remove duplicates
df_out = remove_duplicates(df_all)

# make a xlsx-file, anonymous version
make_cleaned_data(df_out)
