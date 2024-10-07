import os.path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy
import argparse


method_dic = {'mscluster':{'filename':'#Filename','scan':'#Scan','mass':'#ParentMass','rt_time':'#RetTime','cluster':'#ClusterIdx'},
              'falcon':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time','cluster':'cluster'},
              'maracluster':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time','cluster':'cluster'}}

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
    # print("Number of Spectra in mscluster after DB-ID:", len(merged_data))
    return merged_data.groupby('#ClusterIdx').apply(calculate_cluster_purity),merged_data.groupby('#ClusterIdx').size()

def falcon_purity(cluster_results,database_results):
    # Create a merged dataset with left join to include all cluster results and match with database results
    database_results['MzIDFileName'] = database_results['MzIDFileName'].apply(os.path.basename)
    cluster_results['filename'] = cluster_results['filename'].apply(os.path.basename)
    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    # merged_data.to_csv('falcon_merged.csv', index=False)
    largest_cluster_index = merged_data['cluster'].value_counts().idxmax()

    # Print the cluster index of the largest cluster
    # print("Cluster Index of the Largest Cluster:", largest_cluster_index)
    return merged_data.groupby('cluster').apply(calculate_cluster_purity), merged_data.groupby('cluster').size()

def maracluster_purity(cluster_results,database_results):
    database_results['MzIDFileName'] = database_results['MzIDFileName'].apply(os.path.basename)
    cluster_results['filename'] = cluster_results['filename'].apply(os.path.basename)
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

def n10_wrapoer(clustering_results,total_spectra,methods):
    return calculate_n50(clustering_results.groupby(method_dic[methods]['cluster']).size(),total_spectra)

def purity_wapper(clustering_ressults,database_results,methods):
    if methods == 'mscluster':
        return mscluster_purity(clustering_ressults, copy.copy(database_results))
    elif methods == 'falcon':
        return falcon_purity(clustering_ressults, copy.copy(database_results))
    elif methods == 'maracluster':
        return maracluster_purity(clustering_ressults, copy.copy(database_results))

if __name__ == '__main__':
    #read args from command line
    parser = argparse.ArgumentParser(description='Using DB-Search to benchmark')
    parser.add_argument('-c', type=str, required=True, default="cluster_info.tsv", help='input clustering reuslts filename')
    parser.add_argument('-d', type=str, required=True, default="DB_search_demo.tsv",help='input database search results filename')
    parser.add_argument('-t',type=int, required=True,default = 109333, help='Number of MS/MS in the datasets')
    parser.add_argument('-methods', type=str, required=True, default="falcon", help='Clustering methods')
    args = parser.parse_args()
    clustering_filename = args.c
    methods = args.methods
    total_num_spec = args.t
    databse_filename = args.d

    # Load cluster results
    clustering_file_path = '../data/'+methods+'/'+clustering_filename
    database_file_path = '../data/Database-Search/'+databse_filename
    cluster_results = pd.read_csv(clustering_file_path, sep='\t')
    database_results = pd.read_csv(database_file_path, sep='\t')
    database_results['Peptide'].str.replace('I', 'L')

    #calculate n10 number & average purity
    cluster_n10 = n10_wrapoer(cluster_results, total_num_spec,methods)
    print("N10 value for {} clustering results is {}".format(methods,cluster_n10))
    clustering_purity,clustering_size = purity_wapper(cluster_results,database_results,methods)
    avg_purity = sum(clustering_purity) / len(clustering_purity)
    print("Average Purity for {} is {}:".format(methods,avg_purity))
