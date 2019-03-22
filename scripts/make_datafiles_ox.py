import numpy as np
import pandas as pd
from numpy import percentile

datafiles = [
    '../../../Experiments/data/WoT/all_apps_wide_2018-06-19.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-06-20.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-03_reviewed.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-05.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-06.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-10.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-16.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-17.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-19.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-24_reviewed.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-07-25_reviewed.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-08-01_reviewed.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2018-08-09.csv',
    '../../../Experiments/data/WoT/all_apps_wide_2019-03-22.csv',
]

max_error_rate = 10

def make_dataframe(file):
    df = pd.DataFrame(file)

    # take only the relevant columns in data sheet
    df = df[[
            'participant._current_app_name',
            'participant.mturk_worker_id',
            'session.config.selection_method',
            'ox.1.subsession.views',
            'session.code',
            'participant.code',
            'participant.id_in_session',
            'ox.1.player.decision_order',
            'ox.1.player.decision_history',
            'ox.1.player.estimate_kg',
            'ox.1.player.unit',
            'ox.1.player.payoff',
            ]]

    # rename columns
    df.columns = [
        'task', 'mturker', 'method', 'v', 'session', 'code',
        'id', 'order', 'hist', 'guess', 'pound', 'bonus',
    ]

    # drop all rows with nan's
    df = df.dropna()

    return df


def remove_duplicates(df_all):
    duplicates = df_all[df_all['mturker'].duplicated(keep=False)]['mturker'].unique()
    df = df_all.drop_duplicates(['mturker'], keep='first')
    print(len(df_all.guess.values) - len(df.guess.values), 'duplicates removed...\n')
    return df


def error_rate(number, truth):
    return abs(number-truth)/truth


def remove_outliers_error_rate(df_all):
    sessions = df_all.session.unique()
    newdf = pd.DataFrame()
    true_number_of_dots = 1233
    print('\n Now trimming data....')

    # set a lower bound manually
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots

    outliers = []
    above_10 = above_30 = above_50 = above_100 = below = 0
    for session in sessions:
        df_session = df_all[df_all['session'] == session]
        # df_s = df_session[(error_rate(df_session.guess, true_number_of_dots) <= max_error_rate) & (df_session.guess >= lower)]
        # df_out = df_session[(error_rate(df_session.guess, true_number_of_dots) > max_error_rate) | (df_session.guess < lower)]
        # df_s = df_session[(df_session.guess <= upper) & (df_session.guess >= lower)]
        df_s = df_session[(df_session.guess <= upper) & (df_session.guess >= lower)]
        newdf = newdf.append(df_s)

        # stats on outliers:
        df_out = df_session[(df_session.guess > upper) | (df_session.guess < lower)]
        out = df_out.guess.values
        above_10 += len([1 for i in out if error_rate(i,true_number_of_dots) > 10])
        above_30 += len([1 for i in out if error_rate(i,true_number_of_dots) > 30])
        above_50 += len([1 for i in out if error_rate(i,true_number_of_dots) > 50])
        above_100 += len([1 for i in out if error_rate(i,true_number_of_dots) > 100])
        below += len([1 for i in out if i < lower])
        outliers.extend(out)

    print('\noutliers above error rates 10, 30, 50, 100:', above_10, above_30, above_50, above_100)
    print('outliers below truth/10:', below, '\n')
    print(len(df_all.guess.values) - len(newdf.guess.values), 'outliers removed...')
    print(len(newdf.guess.values), 'guesses left ...\nby',
          len(newdf.mturker.unique()), 'unique respondents')
    return newdf


def remove_bad_runs(df):
    # remove experiments which were using the wrong code or made for other
    # purposes
    df = df[df.session != 'j0ssl2ul']  # v=3, history was wrongly initialized
    df = df[df.session != 'gu97xc8e']  # v=1, max was wrongly initialized
    df = df[df.session != 'ka00la5x']  # v=0, history was wrongly initialized
    df = df[df.session != '656wmc92']  # v=3, history was wrongly initialized
    df = df[df.session != '654w33sp']  # v=3, max was wrongly initialized
    df = df[df.session != '8nqnl0a1']  # v=3, too few observations
    print('\nremoving sessions j0ssl2ul, gu97xc8e, ka00la5x, 656wmc92, 8nqnl0a1, and 654w33sp...\n', )
    return df


def make_csvs(df):
    sessions = df.session.unique()
    df_untrimmed_anonymous = pd.DataFrame()
    for session in sessions:
        df_s = df[df['session'] == session]

        # sort with "order"
        df_s = df_s.sort_values('order')

        # make an extra colunm with the true values to emulate the dots data
        df_s['d'] = [1233 for i in range(len(df_s))]

        # drop mturker-column to save anonymous sessions
        df_s = df_s.drop(['mturker'], 1)
        df_s = df_s.drop(['pound'], 1)
        # df_s = df_s.drop(['method'], 1)
        df_untrimmed_anonymous = df_untrimmed_anonymous.append(df_s)

    # reset index
    df_untrimmed_anonymous = df_untrimmed_anonymous.reset_index(drop=True)

    # rearrange position of columns:
    df_untrimmed_anonymous = df_untrimmed_anonymous[['task', 'd', 'method', 'v',
    'session', 'code', 'id', 'order', 'hist', 'guess', 'bonus']]

    # save
    df_untrimmed_anonymous.to_excel('../data/ox.xls')
    print('ox.xls with', len(df_untrimmed_anonymous.guess.values), 'estimates saved.')


# main code
df_all = pd.DataFrame()

for datafile in datafiles:
    df = make_dataframe(pd.read_csv(datafile))
    df_all = df_all.append(df)
print(len(df_all.guess.values), 'guesses found and...\n',
      len(df_all.mturker.unique()), 'unique respondents found...')

# retain only history treatments:
df_all = df_all[df_all.method == 'history']
df_all = df_all[df_all.v != 27]
# remove bad runs and duplicates
df_out = remove_bad_runs(df_all)
# remove duplicates:
df_out = remove_duplicates(df_out)
# make a csv-file for each session, including an anonymous version
make_csvs(df_out)
# trim
remove_outliers_error_rate(df_out)
# remove unconfirmed guesses
# df_out= remove_unconfirmed(df_out)
# make an outlier cleaned csv-file for each session, and an anonymous ditto
# make_trimmed_csvs(df_out)
