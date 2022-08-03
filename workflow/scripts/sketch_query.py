import os, glob
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

df= pd.read_csv(snakemake.input.query_MAG,sep='\t')

path = df.iloc[0]['genome_file']
group = df.iloc[0]['group']


if f"{group}.sketch" in os.listdir("output/reference_sketches"):

    shell(f"bindash sketch --nthreads={snakemake.threads} --kmerlen={snakemake.params.kmer_len} --outfname={snakemake.resources.tmpdir}/{snakemake.wildcards.MAG}.sketch {path} &> {snakemake.log}")
    shell(f"bindash dist --nthreads={snakemake.threads} --outfname={snakemake.output.result} {snakemake.resources.tmpdir}/{snakemake.wildcards.MAG}.sketch output/reference_sketches/{group}.sketch 2>> {snakemake.log}")

else:
    logger.warning("Group was not found in reference! Skipping it . . .")
    shell(f"touch {snakemake.output.result}")