import pandas as pd
import hashlib

raw_datafiles = [
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
            'participant.code',
            'participant.id_in_session',
            'dots.1.player.decision_order',
            'dots.1.player.decision_history',
            'dots.1.player.guess',
            ]]

    # rename columns
    df.columns = [
                  'task',
                  'method',
                  'd',
                  'v',
                  'session',
                  'mturker',
                  'code',
                  'id_in_session',
                  'decision_order',
                  'history',
                  'guess',
    ]
    # drop all rows with nan's
    df = df.dropna()
    return df

def convert_hist_string_to_guess_array(df):
    codes = df.code.values  # these are all the individuals having made a guess
    new_hist_column = []
    for idx, code in enumerate(codes):
        seen_ids = pd.eval(df[df.code == code]['history'].item())
        session = df[df.code == code]['session'].item()
        if seen_ids:
            seen_guesses = [df[(df.id_in_session == g) & (df.session == session)]['guess'].item() for g in seen_ids]
        else:
            seen_guesses = []
        new_hist_column.append(seen_guesses)
    return new_hist_column

def remove_duplicates(df_all):
    images = df_all.d.unique()
    newdf = pd.DataFrame()
    for image in images:
        df_i = df_all[df_all['d'] == image]
        duplicates = df_i[df_i['hashed_turker'].duplicated(keep=False)]['hashed_turker'].unique()
        df_i = df_i.drop_duplicates(['hashed_turker'], keep='first')
        newdf = newdf.append(df_i)
    print(len(df_all.guess.values) - len(newdf.guess.values), 'duplicates removed...')
    return newdf

def make_cleaned_data(df):
    # convert hist-column from string to array of seen guesses:
    df['hist'] = convert_hist_string_to_guess_array(df)

    # make anonymous and drop turker id and code and hist column
    df['hashed_turker'] = hash_simple(df, 'mturker')
    df = df.drop(['mturker', 'code', 'id_in_session', 'history'], 1)

    # remove duplicates:
    df = remove_duplicates(df)

    # sort, rearrange columns, reset index, and save
    df = df.sort_values(by=['method', 'd', 'v', 'session', 'decision_order'])
    df = df[['task', 'method', 'd', 'v', 'session', 'hashed_turker', 'decision_order', 'hist', 'guess']]
    df = df.reset_index(drop=True)
    df.to_excel('../CleanedData/dots.xlsx')
    df.to_csv('../CleanedData/dots.csv')
    print('dots.xlsx and dots.csv saved containing', len(df.guess.values), 'estimates \n \
    from', len(df.hashed_turker.unique()), 'turkers.')


# start here
df_all = pd.DataFrame()
for datafile in raw_datafiles:
    df = make_dataframe(pd.read_csv(datafile))
    df_all = df_all.append(df)

# take only the article-relevant runs:
df_all = df_all[(df_all.method == 'history') | (df_all.method == 'max')]
df_all = df_all[df_all.v != 27]
print(len(df_all.guess.values), 'estimates found and...\n',
      len(df_all.mturker.unique()), 'unique respondents found...')

# make a xlsx-file, anonymous version
make_cleaned_data(df_all)
