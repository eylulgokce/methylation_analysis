import pandas as pd
import numpy as np
import os
import gzip
import time
import sys
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import entropy, ks_2samp
from scipy.special import rel_entr
from scipy.spatial.distance import jensenshannon


############ Divergence ##############

def kl_divergence(p, q):
    return np.sum(rel_entr(p, q))

def jensen_shannon_divergence(p, q):
    m = 0.5 * (p + q)
    return 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)

def geometric_jsd(p, q):
    jsd = jensen_shannon_divergence(p, q)
    return np.sqrt(jsd)

def ks_test(p, q):
    ks_result = ks_2samp(p, q)
    return ks_result.statistic, ks_result.pvalue

############# Smoothing ###################

def laplace_smoothing(data, alpha=1e-10):
    return (data + alpha) / (np.sum(data) + alpha * len(data))

#################### ANALYZE ########################

def process_and_analyze(file_path, output_path):
    try:
        with gzip.open(file_path, 'rt') as f:
            df = pd.read_csv(f, sep='\t', header=None, names=['chr', 'start', 'end', 'percentage', 'methylated', 'unmethylated'])

        results_list = []

        for index, row in df.iterrows():
            start_time = time.time()

            # Prepare data
            data = np.array([row['methylated'], row['unmethylated']]) 
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

        # check if the output directory exists
        #os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        results_df.to_csv(output_path, index=False, compression='gzip')

        return results_df
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None


########################### MAIN ###############################

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("ARGS is not equal to 3! USE: python run.py <input_file.gz> <output_file.csv.gz>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_and_analyze(input_file, output_file)
