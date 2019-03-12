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
    '../../../Experiments/wot/data/ox/all_apps_wide_2019-03-08.csv',
    '../../../Experiments/wot/data/ox/all_apps_wide_2019-03-11.csv',
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
        'task', 'mturker', 'method', 'views', 'session', 'code',
        'id', 'order', 'hist', 'guess', 'pound', 'bonus',
    ]

    # drop all rows with nan's
    df = df.dropna()

    return df


def remove_duplicates(df_all):
    duplicates = df_all[df_all['mturker'].duplicated(keep=False)]['mturker'].unique()
    df = df_all.drop_duplicates(['mturker'], keep='first')
    print(len(df_all.guess.values) - len(df.guess.values), 'duplicates removed...')
    print(len(df.guess.values), 'guesses left ...\nand',
          len(df.mturker.unique()), 'unique respondents left...')
    print('containing', len(df[df['method'] == 'history'].guess.values),
          'subjects in the history condition, ...\n',
          len(df[df['method'] == 'max'].guess.values),
          'subjects in the max condition, ...\n',
          len(df[df['method'] == 'min'].guess.values),
          'subjects in the min condition...\n')
    return df


def error_rate(number, truth):
    return abs(number-truth)/truth


def remove_outliers_error_rate(df_all):
    sessions = df_all.session.unique()
    newdf = pd.DataFrame()
    true_number_of_dots = 1233
    print('\n Now trimming data....')

    # set a lower bound manually (get rid of the cranks)
    lower = true_number_of_dots/10
    upper = max_error_rate*true_number_of_dots + true_number_of_dots
    print('upper bound =', upper)
    print('lower bound =', lower)

    outliers = []
    above_10 = above_30 = above_50 = above_100 = 0
    for session in sessions:
        df_session = df_all[df_all['session'] == session]
        # df_s = df_session[(error_rate(df_session.guess, true_number_of_dots) <= max_error_rate) & (df_session.guess >= lower)]
        # df_out = df_session[(error_rate(df_session.guess, true_number_of_dots) > max_error_rate) | (df_session.guess < lower)]
        df_s = df_session[(df_session.guess <= upper) & (df_session.guess >= lower)]
        df_out = df_session[(df_session.guess > upper) | (df_session.guess < lower)]

        out = df_out['guess'].values
        above_10 += len([1 for i in out if error_rate(i,true_number_of_dots) > 10])
        above_30 += len([1 for i in out if error_rate(i,true_number_of_dots) > 30])
        above_50 += len([1 for i in out if error_rate(i,true_number_of_dots) > 50])
        above_100 += len([1 for i in out if error_rate(i,true_number_of_dots) > 100])
        outliers.extend(out)
        newdf = newdf.append(df_s)

    print(outliers, 'number of outliers:', len(outliers))
    print('\noutliers above error rates 10, 30, 50, 100:', above_10, above_30, above_50, above_100)
    print(len(df_all.guess.values) - len(newdf.guess.values), 'outliers removed...')
    print(len(newdf.guess.values), 'guesses left ...\nand',
          len(newdf.mturker.unique()), 'unique respondents left...')
    print('containing', len(newdf[newdf['method'] == 'history'].guess.values),
          'subjects in the history condition, ...\n',
          len(newdf[newdf['method'] == 'max'].guess.values),
          'subjects in the max condition, ...\n',
          len(newdf[newdf['method'] == 'min'].guess.values),
          'subjects in the min condition...\n')
    return newdf


def remove_outliers(df_all):
    sessions = df_all.session.unique()
    newdf = pd.DataFrame()

    above_10 = above_30 = above_50 = above_70 = 0
    for session in sessions:
        df_s = df_all[df_all['session'] == session]
        views = df_s['views'].unique().item()
        method = df_s['method'].unique().item()
        print('\nSession', session, views, method)
        data = df_s['guess'].values

        # calculate interdecile range
        q10, q90 = percentile(data, 10), percentile(data, 90)
        ridr = (q90 - q10)/np.median(data)
        idr = q90 - q10
        print('deciles: 10th=%.3f, 90th=%.3f, IDR=%.3f, RIDR=%.3f, median=%.3f' % (q10, q90, idr, ridr,np.median(data)))

        # calculate the outlier cutoff
        cut_off = idr * 3
        upper = q90 + cut_off

        # set a lower bound manually (get rid of the cranks)
        lower = 1233/10
        print('lower bound =', lower, 'upper bound=', upper)

        # identify outliers
        outliers = [x for x in data if x < lower or x > upper]
        print('Identified outliers: %d' % len(outliers), 'out of', len(data))
        print(outliers)
        above_10 += len([1 for i in outliers if error_rate(i,1233) > 10])
        above_30 += len([1 for i in outliers if error_rate(i,1233) > 30])
        above_50 += len([1 for i in outliers if error_rate(i,1233) > 50])
        above_70 += len([1 for i in outliers if error_rate(i,1233) > 70])

        # remove outliers
        df_s = df_s[(df_s.guess >= lower) & (df_s.guess <= upper)]
        newdf = newdf.append(df_s)

    print('outliers above error rates 10, 30, 50, 70:', above_10, above_30, above_50, above_70)
    print(len(df_all.guess.values) - len(newdf.guess.values), 'outliers removed...')
    print(len(newdf.guess.values), 'guesses and',
          len(newdf.mturker.unique()), 'unique respondents left...')
    print('containing', len(newdf[newdf['method'] == 'history'].guess.values),
          'subjects in the history condition, ...\n',
          len(newdf[newdf['method'] == 'max'].guess.values),
          'subjects in the max condition, ...\n',
          len(newdf[newdf['method'] == 'min'].guess.values),
          'subjects in the min condition...\n')
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
    print(len(df.guess.values), 'guesses left ...\nand ',
          len(df.mturker.unique()), 'unique respondents left...')
    return df


def make_csvs(df):
    # make a txt file showing session names
    sessions = df.session.unique()
    # np.savetxt("../published/data/ox/sessions.csv", sessions,
    #            fmt='%s', delimiter="")
    print(len(sessions), 'session names saved in sessions.csv...')

    # save a new csv-file for each session and one for all
    df_untrimmed = pd.DataFrame()
    df_untrimmed_anonymous = pd.DataFrame()
    for session in sessions:
        df_s = df[df['session'] == session]

        # sort with "order"
        df_s = df_s.sort_values('order')

        # make an extra colunm with the true values to emulate the dots data
        df_s['dots'] = [1233 for i in range(len(df_s))]

        # save the individual sessions
        # df_s.to_csv('../published/data/ox/' + str(session) + '_untrimmed.csv')
        # print(str(session)+'_untrimmed.csv with', len(df_s.guess.values), 'guesses saved')
        df_untrimmed = df_untrimmed.append(df_s)

        # drop mturker-column to save anonymous sessions
        df_s = df_s.drop(['mturker'], 1)
        df_s = df_s.drop(['method'], 1)
        # df_s.to_csv('../published/data/ox/' + str(session) + '_untrimmed_anonymous.csv')
        # print(str(session)+'_untrimmed_anonymous.csv with', len(df_untrimmed.guess.values), 'guesses saved')
        df_untrimmed_anonymous = df_untrimmed_anonymous.append(df_s)

    print(len(sessions), 'threads saved as _untrimmed.csv\'s and _untrimmed_anonymous.csv\'s')

    # reset index
    df_untrimmed = df_untrimmed.reset_index(drop=True)
    df_untrimmed_anonymous = df_untrimmed_anonymous.reset_index(drop=True)

    # save in one untrimmed and one untrimmed_anonymous csv-file
    # df_untrimmed.to_csv('../published/data/ox/all_ox_untrimmed.csv')
    # print('all_ox_untrimmed.csv with', len(df_untrimmed.guess.values), 'raw guesses saved...')
    # df_untrimmed_anonymous.to_csv('../published/data/ox/all_ox_untrimmed_anonymous.csv')
    df_untrimmed_anonymous.to_excel('../data/ox/ox.xlsx')
    print('all_ox_untrimmed_anonymous.csv with', len(df_untrimmed_anonymous.guess.values), 'guesses saved...')
    print('containing', len(df_untrimmed[df_untrimmed['method'] == 'history'].guess.values),
          'subjects in the history condition, ...\n',
          len(df_untrimmed[df_untrimmed['method'] == 'max'].guess.values),
          'subjects in the max condition, ...\n',
          len(df_untrimmed[df_untrimmed['method'] == 'min'].guess.values),
          'subjects in the min condition...\n')


def make_trimmed_csvs(df):
    sessions = df.session.unique()

    # save a new csv-file for each session and one for all
    df_trimmed = pd.DataFrame()
    df_trimmed_anonymous = pd.DataFrame()
    for session in sessions:
        df_s = df[df['session'] == session]

        # sort with "order"
        df_s = df_s.sort_values('order')

        # make an extra colunm with the true values to emulate the dots data
        df_s['dots'] = [1233 for i in range(len(df_s))]

        # save the outlier-trimmed data
        # df_s.to_csv('../published/data/ox/' + str(session) + '_trimmed.csv')
        # print(str(session)+'_trimmed.csv with', len(df_s.guess.values), 'guesses saved')
        df_trimmed = df_trimmed.append(df_s)

        # drop mturker-column to make data anonymous and save
        df_s = df_s.drop(['mturker'], 1)
        # df_s.to_csv('../published/data/ox/' + str(session) + '_trimmed_anonymous.csv')
        # print(str(session)+'_trimmed_anonymous.csv with', len(df_s.guess.values), 'guesses saved')
        df_trimmed_anonymous = df_trimmed_anonymous.append(df_s)

    print(len(sessions), 'threads saved as _trimmed.csv\'s and _trimmed_anonymous.csv\'s')

    # reset index
    df_trimmed = df_trimmed.reset_index(drop=True)
    df_trimmed_anonymous = df_trimmed_anonymous.reset_index(drop=True)

    # save in one clean and one anonymous csv-file
    # df_trimmed.to_csv('../published/data/ox/all_ox_trimmed.csv')
    print('all_ox_trimmed.csv with', len(df_trimmed.guess.values), 'guesses saved')
    # df_trimmed_anonymous.to_csv('../published/data/ox/all_ox_trimmed_anonymous.csv')
    print('all_ox_trimmed_anonymous.csv with', len(df_trimmed_anonymous.guess.values), 'guesses saved')
    print('containing', len(df_trimmed[df_trimmed['method'] == 'history'].guess.values),
          'subjects in the history condition, ...\n',
          len(df_trimmed[df_trimmed['method'] == 'max'].guess.values),
          'subjects in the max condition, ...\n',
          len(df_trimmed[df_trimmed['method'] == 'min'].guess.values),
          'subjects in the min condition...\n')

# main code
df_all = pd.DataFrame()

for datafile in datafiles:
    df = make_dataframe(pd.read_csv(datafile))
    df_all = df_all.append(df)
print(len(df_all.guess.values), 'guesses found and...\n',
      len(df_all.mturker.unique()), 'unique respondents found...')

# retain only history treatments:
df_all = df_all[df_all.method == 'history']
# remove bad runs and duplicates
df_out = remove_bad_runs(df_all)
# remove duplicates:
df_out = remove_duplicates(df_out)
# make a csv-file for each session, including an anonymous version
make_csvs(df_out)
# trim
df_out = remove_outliers_error_rate(df_out)
# remove unconfirmed guesses
# df_out= remove_unconfirmed(df_out)
# make an outlier cleaned csv-file for each session, and an anonymous ditto
make_trimmed_csvs(df_out)
