import os.path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy
import argparse

def calculate_cluster_purity(cluster):
    total_spectra = len(cluster)
    if total_spectra == 1:
        return 1

    identified_peptides_counts = cluster.groupby('Peptide').size()
    max_count = identified_peptides_counts.max() if not identified_peptides_counts.empty else 0

    max_fraction = max_count/total_spectra

    return max_fraction

def mscluster_purity(cluster_results,database_results):
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzML')
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
    # Create a merged dataset with left join to include all cluster results and match with database results
    merged_data = cluster_results.merge(database_results, left_on=['#Filename', '#Scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    merged_data = merged_data[merged_data['Peptide'] != 'unknown']
    print("Number of Spectra in mscluster after DB-ID:", len(merged_data))
    return merged_data.groupby('#ClusterIdx').apply(calculate_cluster_purity),merged_data.groupby('#ClusterIdx').size()

def falcon_purity(cluster_results,database_results):
    # Create a merged dataset with left join to include all cluster results and match with database results
    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    merged_data.to_csv('falcon_merged.csv', index=False)
    largest_cluster_index = merged_data['cluster'].value_counts().idxmax()

    # Print the cluster index of the largest cluster
    print("Cluster Index of the Largest Cluster:", largest_cluster_index)
    return merged_data.groupby('cluster').apply(calculate_cluster_purity), merged_data.groupby('cluster').size()

def maracluster_purity(cluster_results,database_results):
    database_results['MzIDFileName'] = database_results['MzIDFileName'].apply(os.path.basename)
    cluster_results['filename'] = cluster_results['filename']+'.mzML'
    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    return merged_data.groupby('cluster').apply(calculate_cluster_purity), merged_data.groupby('cluster').size()

def assign_size_range_bin(size):
    if size == 1:
        return '1'
    else:
        upper_limit = 2
        while size > upper_limit:
            upper_limit *= 2
        return f"{upper_limit//2+1} to {upper_limit}"

def online_falcon_purity(cluster_results,database_results):
    # Create a merged dataset with left join to include all cluster results and match with database results
    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    merged_data.to_csv('online_merged.csv', index=False)
    largest_cluster_index = merged_data['cluster'].value_counts().idxmax()

    # Print the cluster index of the largest cluster
    print("Cluster Index of the Largest Cluster:", largest_cluster_index)
    return merged_data.groupby('cluster').apply(calculate_cluster_purity), merged_data.groupby('cluster').size()

def calculate_weighted_average_purity(purity, cluster_size):
    if cluster_size.sum() == 0:
        return 0
    weighted_purity = (purity * cluster_size).sum() / cluster_size.sum()
    return weighted_purity

def calculate_n50(cluster_size, total_spectra):
    sorted_sizes = np.sort(cluster_size)[::-1]  # Sort cluster sizes in descending order
    #total_spectra = sorted_sizes.sum()
    cumulative_sum = 0
    n50 = None
    for size in sorted_sizes:
        cumulative_sum += size
        if cumulative_sum >= total_spectra *0.1:
            n50 = size
            break
    return n50

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Using ClassyFIre results to benchmark the network')
    parser.add_argument('-c', type=str, required=True, default="cluster_info.tsv", help='input clustering reuslts filename')
    parser.add_argument('-d', type=str, required=True, default="demo.tsv",help='input database search results filename')
    parser.add_argument('-t',type=int, required=True,default = None, help='Number of MS/MS in the datasets')
    parser.add_argument('-methods', type=str, required=True, default="Falcon", help='Clustering methods')
    args = parser.parse_args()
    # Load cluster results
    maracluster_results= pd.read_csv('../data/PXD023047/maracluster/MaRaCluster_processed.clusters_p2_enriched.tsv', sep='\t')
    falcon_results = pd.read_csv('/home/user/LabData/XianghuData/Falcon_Cluster_Benchmark/online_clustering_new_0.25/output_summary/cluster_info.tsv',sep='\t')
    online_falcon_results = pd.read_csv('/home/user/research/online_clustering/work_0.25/cluster_info.tsv',sep='\t')
    # online_falcon_results = pd.read_csv('../data/online_clustering/online_cluster_into.tsv',sep='\t')
    online_falcon_results = online_falcon_results[online_falcon_results['cluster'] != -1]
    database_results = pd.read_csv('../data/online_clustering/online_new_filtered.tsv', sep='\t')  # Adjust file path and format accordingly
    #print("oroginal number of Spectra in mscluster:",len(mscluster_results))
    # mscluster_n50 = calculate_n50(mscluster_results.groupby('#ClusterIdx').size(),109333)
    falcon_n50 = calculate_n50(falcon_results.groupby('cluster').size(),6909884)
    online_falcon_n50 = calculate_n50(online_falcon_results.groupby('cluster').size(),6909884)
    #maracluster_n50 = calculate_n50(maracluster_results.groupby('cluster').size(),286410)


    # print("N50 value for MSCluster:", mscluster_n50)
    print("N50 value for Falcon:", falcon_n50)
    print("N50 value for Online Falcon:", online_falcon_n50)
    #print("N50 value for MaRaCluster:", maracluster_n50)

    # mscluster_results = pd.read_csv('../data/results/nf_output/clustering/clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    # falcon_results = pd.read_csv('../data/cluster_info.tsv', sep='\t')  # Adjust file path and format accordingly
    # online_falcon_results = pd.read_csv('../data/online_cluster_info.tsv', sep='\t')  # Adjust file path and format accordingly
    # maracluster_results= pd.read_csv('./processed_clusters.tsv', sep='\t')
    # database_results = pd.read_csv('./filtered.tsv', sep='\t')



    #mscluster_purity,mscluster_size = mscluster_purity(mscluster_results,copy.copy(database_results))
    falcon_purity,falcon_size =falcon_purity(falcon_results,copy.copy(database_results))
    online_falcon_purity,online_falcon_size =online_falcon_purity(online_falcon_results,copy.copy(database_results))
    #maracluster_purity,maracluster_size = maracluster_purity(maracluster_results,copy.copy(database_results))
    database_results['Peptide'].str.replace('I', 'L')


    # mscluster_weighted_avg_purity = calculate_weighted_average_purity(mscluster_purity, mscluster_size)
    #falcon_weighted_avg_purity = calculate_weighted_average_purity(falcon_purity, falcon_size)
    falcon_weighted_avg_purity = sum(falcon_purity)/len(falcon_purity)
    #online_falcon_weighted_avg_purity = calculate_weighted_average_purity(online_falcon_purity, online_falcon_size)
    online_falcon_weighted_avg_purity =sum(online_falcon_purity)/len(online_falcon_purity)
    #maracluster_weighted_avg_purity = calculate_weighted_average_purity(maracluster_purity, maracluster_size)


    # print("Weighted Average Purity for MSCluster:", mscluster_weighted_avg_purity)
    print("Weighted Average Purity for Falcon:", falcon_weighted_avg_purity)
    print("Weighted Average Purity for Online Falcon:", online_falcon_weighted_avg_purity)
