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
os.makedirs(snakemake.output.batch_dir, exist_ok=True)

paths = snakemake.config["MAGS_path_table"]

paths_df = pd.read_csv(paths, sep="\t")[['path','group']]


# Split based on batch_size in random batches

groups = paths_df.group.unique()
for group in groups:
    family_df = paths_df.loc[paths_df["group"] == group]
    names = (
        family_df.path.str.rsplit("/", 1, expand=True)[1]
        .str.rsplit(".", 2, expand=True)[0]
        .values.tolist()
    )
    output_df = pd.DataFrame()
    output_df["F"] = family_df["path"].values.tolist()
    group = group.lstrip("f__")
    output_df.to_csv(f"{snakemake.output.batch_dir}/{group}.tsv", sep='\t',index=False,header=False)
    logging.info(f"Created batch file {snakemake.output.batch_dir}/{group}.tsv")