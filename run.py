import pandas as pd
import numpy as np
import os
import gzip
import time
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import entropy, ks_2samp
from scipy.special import rel_entr
from scipy.spatial.distance import jensenshannon

# Function to calculate KL divergence
def kl_divergence(p, q):
    return np.sum(rel_entr(p, q))

def geometric_jsd(p, q):
    jsd = jensenshannon(p, q)
    return np.sqrt(jsd)


def ks_test(p, q):
    ks_result = ks_2samp(p, q)
    return ks_result.statistic, ks_result.pvalue

def laplace_smoothing(data, alpha=1e-10):
    return (data + alpha) / (np.sum(data) + alpha * len(data))

def process_and_analyze(file_path, output_base_dir):
    with gzip.open(file_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None, names=['chr', 'start', 'end', 'percentage', 'methylated', 'unmethylated'])

    results_list = []

    for index, row in df.iterrows():
        start_time = time.time()

        # Prepare data
        data = np.array([row['methylated'], row['unmethylated']]) # Laplace increases both methylated and unmethylated counts
        smoothed_data = laplace_smoothing(data)
        p = np.array([smoothed_data[0], smoothed_data[1]])
        q = np.array([smoothed_data[1], smoothed_data[0]])

        # Normalize
        p /= p.sum()
        q /= q.sum()

        # Calculate divergence
        ent = entropy(p)
        kl = kl_divergence(p, q)
        js = jensenshannon(p, q)
        gjs = geometric_jsd(p, q)
        ks_stat, ks_pvalue = ks_test(p, q)
        end_time = time.time()
        time_elapsed = end_time - start_time

        results_list.append([row['chr'], row['start'], row['end'], row['percentage'], row['methylated'], row['unmethylated'], ent, kl, js, gjs, ks_stat, ks_pvalue, time_elapsed])

    results_df = pd.DataFrame(results_list, columns=['chr', 'start', 'end', 'percentage', 'methylated', 'unmethylated', 'entropy', 'relative_entropy', 'jsd', 'geometric_jsd', 'ks_stat', 'ks_pvalue', 'time'])

    sample_dir = os.path.join(output_base_dir, os.path.basename(os.path.dirname(file_path)))
    os.makedirs(sample_dir, exist_ok=True)

    output_filename = os.path.join(sample_dir, f"output_{os.path.basename(file_path).replace('.gz', '')}.csv.gz")
    results_df.to_csv(output_filename, index=False, compression='gzip')

    return results_df

base_dir = '/shares/grossniklaus.botinst.uzh/eharputluoglu/datasets_gz'
output_base_dir = '/shares/grossniklaus.botinst.uzh/eharputluoglu/output'

files_to_process = []
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith('.cov.gz'):
            files_to_process.append(os.path.join(root, file))

# Process files in parallel
with ProcessPoolExecutor(max_workers=10) as executor:
    futures = [executor.submit(process_and_analyze, file_path, output_base_dir) for file_path in files_to_process]
    for future in futures:
        future.result()
