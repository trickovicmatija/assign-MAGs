import os, glob
from turtle import distance
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool


from snakemake.shell import shell
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

# Start script

logger.info('Starting preparing input data!')

distance_df = pd.read_csv(snakemake.input.distance,sep='\t')

distance_df['hash_found'] = distance_df['hash_ratio'].str.split("/",expand=True)[0].astype('int')
distance_df['hash_total'] = distance_df['hash_ratio'].str.split("/",expand=True)[1].astype('int')
distance_df['hash_ratio'] = distance_df['hash_found']/distance_df['hash_total']

group_def = pd.read_csv(snakemake.config['group_def'],sep='\t')
group_def.subspecies = group_def.subspecies.astype('str').str.rjust(10,'0')

distance_df['target_genome'] = distance_df['target'].str.rsplit("/",1,expand=True)[1].str.rsplit(".",2,expand=True)[0]

distance_df = distance_df.merge(group_def,left_on='target_genome',right_on='genome',how='left')
logger.info('Finished preparing input data!')

output_df = pd.DataFrame()
short_df = pd.DataFrame()

def assign_group(query_MAG):
    short_output_dict = defaultdict()
    temp_df = distance_df.loc[distance_df['query'] == query_MAG]  
    subsp_dist = defaultdict()
    for subsp in temp_df.subspecies.unique():
        df = temp_df.loc[temp_df['subspecies'] == subsp]
        df = df.loc[df['hash_found'] > int(snakemake.config['min_hash'])]
        if df.shape[0] > 0:
            subsp_dist[subsp] = df['dist'].median()
    dist_df = pd.DataFrame.from_dict(subsp_dist,orient='index')
    dist_df.columns = ['median_dist']
    dist_df = dist_df.reset_index()
    dist_df = dist_df.rename(columns={'index':'subsp'})
    dist_df['query'] = query_MAG

    chosen_subsp = dist_df.sort_values(by='median_dist',ascending=True).iloc[0]['subsp']
    short_output_dict[query_MAG] = chosen_subsp
    short_output_df = pd.DataFrame.from_dict(short_output_dict,orient='index')
    short_output_df.columns = ['subsp']
    short_output_df = short_output_df.reset_index()
    short_output_df = short_output_df.rename(columns={'index':'query'})
    short_output_df['median_dist'] = dist_df.loc[dist_df['subsp'] == chosen_subsp]['median_dist'].iloc[0]

    logger.info(f"Finished assigning group for {query_MAG}")

    return dist_df, short_output_df

all_queries = distance_df['query'].unique().tolist()
logger.info('Starting assigning groups!')

for query in all_queries:
    dist_df, short_output_df = assign_group(query)
    output_df = pd.concat([output_df,dist_df])
    short_df = pd.concat([short_df, short_output_df])

cols = list(output_df.columns)
cols = [cols[-1]] + cols[:-1]
output_df = output_df[cols]

output_df.to_csv(snakemake.output.full_output,sep='\t',index=False)

short_df.to_csv(snakemake.output.short_output, sep='\t', index=False)

logger.info('Finished assigning groups!')