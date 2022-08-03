import os, glob
import pandas as pd
from collections import defaultdict

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
output_df = pd.DataFrame()

for distance_df in snakemake.input.all_input:
    MAG = distance_df.split("/")[-1].split(".")[0]
    if os.path.getsize(distance_df) > 0:
        df = pd.read_csv(distance_df,sep='\t',names=['query','target','dist','p_value','hash_ratio'],header=None)
        output_df = pd.concat([output_df,df],axis=0)
    else:
        logger.warning(f"Table for MAG {MAG} is empty!")

output_df.to_csv(snakemake.output.distance_df,sep='\t',index=False)