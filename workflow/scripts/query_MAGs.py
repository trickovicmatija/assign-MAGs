import os
import pandas as pd

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
os.makedirs(snakemake.output.query_dir, exist_ok=True)

df_path = snakemake.config['query_paths']
df = pd.read_csv(df_path,sep='\t')

for MAG in df.Genome:
    temp_df = df.query("Genome == @MAG")
    path = temp_df['genome_file'].iloc[0]
    if os.path.exists(path):
        temp_df.to_csv(f"{snakemake.output.query_dir}/{MAG}.tsv", sep='\t',index=False)
    else:
        logger.warning(f"Path not found for {path}")
