import pysam
import pandas as pd
import os
from glob import glob

def load_predictions(prediction_dir):
    readid_to_barcode = {}
    for file in glob(os.path.join(prediction_dir, "barcode_predictions_*")):
        df = pd.read_csv(file)
        for _, row in df.iterrows():
            readid_to_barcode[row["read_id"]] = row["predicted_barcode"]
    return readid_to_barcode

def split_bam_by_barcode(bam_path, predictions, output_prefix):
    barcode_to_writer = {}
    bam_in = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    
    def get_writer(barcode):
        if barcode not in barcode_to_writer:
            out_path = f"{output_prefix}_barcode{barcode}.bam"
            barcode_to_writer[barcode] = pysam.AlignmentFile(out_path, "wb", header=bam_in.header)
        return barcode_to_writer[barcode]
    
    for read in bam_in.fetch(until_eof=True):
        read_id = read.get_tag("pi") if read.has_tag("pi") else read.query_name
        barcode = predictions.get(read_id, -1)
        writer = get_writer(barcode)
        writer.write(read)

    for writer in barcode_to_writer.values():
        writer.close()

    bam_in.close()

# Usage
prediction_dir = "./predictions"
bam_path = "./mapped.bam"
output_prefix = "prefix"
predictions = load_predictions(prediction_dir)
split_bam_by_barcode(bam_path, predictions, output_prefix)
