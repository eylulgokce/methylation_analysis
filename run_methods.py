import pandas as pd
import numpy as np
import glob
import os
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from scipy.stats import entropy, ks_2samp, f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='methylation_analysis_3.log', 
    filemode='w'                        
)

def calculate_divergences_and_ks(df, smoothing=1e-5):
    """Calculates normalized divergences (JSD, KL, GJS, SGJS) and KS test for methylation data.

    Args:
        df (pandas.DataFrame): DataFrame with 'methylated' and 'unmethylated' columns.
        smoothing (float, optional): Smoothing factor to add to counts. Defaults to 1e-10.

    Returns:
        tuple: Normalized divergences (JSD, KL, GJS, SGJS) and KS test statistic and p-value.
    """

    # Filter out rows with zero coverage
    #df = df[df['coverage'] > 0].copy()

    # Drop rows where both methylated and unmethylated are zero 
    df = df[(df['methylated'] != 0) | (df['unmethylated'] != 0)].copy()

    if len(df) == 0:
        return None, None, None, None, None, None
    
    df["methylated"] += smoothing
    df["unmethylated"] += smoothing
    
    # Calculate probabilities
    total = df["methylated"] + df["unmethylated"]
    p = df["methylated"] / total
    q = df["unmethylated"] / total
    
    # Divergence calculations (with normalizations)
    normalized_js_divergence = float(jensenshannon(p, q, base=2))  # JSD

    normalized_kl_divergence = float(entropy(p, q, base=2) - entropy(p, base=2))  # KLD

    # Geometric Jensen-Shannon Divergence
    normalized_gjs_divergence = float(np.sqrt(
        0.5 * (np.sqrt(entropy(p, q, base=2)) + np.sqrt(entropy(q, p, base=2)))
    ) / np.sqrt(np.log(2)))
    
    # Symmetric Geometric Jensen-Shannon Divergence
    m = (p + q) / 2
    normalized_sgjs_divergence = float( 0.5 * (
        np.sqrt(jensenshannon(p, m, base=2)) + np.sqrt(jensenshannon(q, m, base=2))
    ) / np.sqrt(np.log(2)))

    # Kolmogorov-Smirnov Test
    ks_statistic, ks_pvalue = ks_2samp(df["methylated"], df["unmethylated"])  # KST

    return normalized_js_divergence, normalized_kl_divergence, normalized_gjs_divergence, normalized_sgjs_divergence, ks_statistic, ks_pvalue

def read_bismark_file(filename):
    """Reads Bismark methylation data from a file."""
    column_names = ["chr", "start", "end", "coverage", "methylated", "unmethylated"]
    df = pd.read_csv(filename, sep='\t', header=None, names=column_names, compression='gzip')

    # Convert columns to correct data types
    df['start'] = pd.to_numeric(df['start'])
    df['end'] = pd.to_numeric(df['end'])
    df['coverage'] = pd.to_numeric(df['coverage'])
    df['methylated'] = pd.to_numeric(df['methylated'])
    df['unmethylated'] = pd.to_numeric(df['unmethylated'])

    return df

def clean_data(results_df):
    """Cleans the results DataFrame by removing file extensions from sample names."""
    results_df['Sample'] = results_df['Sample'].astype(str).str.replace('.bedgraph.gz', '', regex=False)
    return results_df


# Data processing pipeline
data_directory = "/shares/grossniklaus.botinst.uzh/dkt/projects/meth1000/analysis/09_split_cov_chr/output"
results = []

# Get all sample subdirectories
sample_directories = [d for d in os.listdir(data_directory) if os.path.isdir(os.path.join(data_directory, d))]

# Iterate through files
for sample_dir in sample_directories:
    sample_path = os.path.join(data_directory, sample_dir)

    logging.info(f"Processing sample directory: {sample_path}")

    # Iterate through files in each sample directory
    for file_path in glob.glob(os.path.join(sample_path, "*_chr_*.cov.gz")):
        try:
            filename = os.path.basename(file_path)
            parts = filename.split("_")

            context_type = parts[0]
            sample_name = parts[1] + "_" + parts[2].split(".")[0] 
            chromosome = parts[-1].split(".")[0]

            #logging.info(f"Processing file: {filename}")



            df = read_bismark_file(file_path)
            divergences = calculate_divergences_and_ks(df)
    
            results.append(
                {
                    "Sample": sample_name,
                    "Chromosome": chromosome,
                    "Context": context_type,
                    "JS Divergence": divergences[0],  
                    "KL Divergence": divergences[1],
                    "GJS Divergence": divergences[2],
                    "SGJS Divergence": divergences[3],
                    "KS Statistic": divergences[4],
                    "KS P-value": divergences[5]
                }
            )
    
        except Exception as e:
            logging.error(f"Error processing file {filename}: {e}")

# Create a DataFrame from results
results_df = pd.DataFrame(results)

# Clean data
#results_df = clean_data(results_df)

# Save the results to a CSV file
results_df.to_csv("/shares/grossniklaus.botinst.uzh/eharputluoglu/calculation_output/methylation_analysis_results_Smooting_1e-5.csv", index=False)