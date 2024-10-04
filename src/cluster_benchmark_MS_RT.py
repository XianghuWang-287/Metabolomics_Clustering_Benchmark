import os
from itertools import combinations
from pyteomics import mzml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pymzml
from  xml.etree.ElementTree import ParseError
import networkx as nx
from collections import Counter
from matplotlib.colors import LogNorm
import pickle
import copy
from tqdm import tqdm
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
import argparse

method_dic = {'mscluster':{'filename':'#Filename','scan':'#Scan','mass':'#ParentMass','rt_time':'#RetTime','cluster':'#ClusterIdx'},
              'falcon':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time','cluster':'cluster'},
              'maracluster':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time','cluster':'cluster'}}

def get_ms1_data(mzml_file):
    ms1_data = []

    run = pymzml.run.Reader(mzml_file, build_index_from_scratch=False)



    try:
        for idx, spectrum in enumerate(run):
            if spectrum.ms_level == 2:
                precursors = spectrum.selected_precursors
                if len(precursors) != 1:
                    exception = NotImplementedError(
                        "Expected one precursors but got {} while checking index {}.".format(len(precursors),
                                                                                             spectrum.index))
                    print(exception)
                    print(mzml_file)
                    continue
                rt = float(spectrum.scan_time_in_minutes()) * 60
                mass =precursors[0]['mz']
                scan_number = int(spectrum['id'])
                ms1_data.append((mass, rt, scan_number))
    except Exception as e:
            print(e)

    return ms1_data

def are_matching_spectra(filename, scan1, scan2,matching_pairs_set):
    return (filename, scan1, scan2) in matching_pairs_set

def optimized_create_matching_network(cluster, method):
    G = nx.Graph()

    # Precompute the node names and add them to the graph
    node_attrs = {
        f"{row[method_dic[method]['filename']]}_{row[method_dic[method]['scan']]}": {"filename": row[method_dic[method]['filename']]}
        for _, row in cluster.iterrows()
    }
    G.add_nodes_from(node_attrs.items())

    # Precompute mass and rt_time for efficient access
    specs = cluster.apply(lambda row: (
        f"{row[method_dic[method]['filename']]}_{row[method_dic[method]['scan']]}",
        row[method_dic[method]['mass']],
        row[method_dic[method]['rt_time']]
    ), axis=1).tolist()

    # Create edges based on conditions
    edges = [
        (spec1[0], spec2[0]) for spec1, spec2 in combinations(specs, 2)
        if spec1[0].split('_')[0] == spec2[0].split('_')[0]  # Ensure filenames are the same
           and abs(spec1[1] - spec2[1]) <= 0.01 and abs(spec1[2] - spec2[2]) <= 0.5
    ]
    G.add_edges_from(edges)

    return G
def create_matching_network(cluster, method):
    G = nx.Graph()

    row_indices = range(len(cluster))
    row_combinations = combinations(row_indices, 2)
    for _, spec in cluster.iterrows():
        node_a = f"{spec[method_dic[method]['filename']]}_{spec[method_dic[method]['scan']]}"
        G.add_node(node_a, filename=spec[method_dic[method]['filename']])
    for combo in row_combinations:
        index1, index2 = combo
        spec1 = cluster.iloc[index1]
        spec2 = cluster.iloc[index2]
        if spec1[method_dic[method]['filename']] == spec2[method_dic[method]['filename']]:
            if abs(spec1[method_dic[method]['mass']]-spec2[method_dic[method]['mass']])<=0.01 and abs(spec1[method_dic[method]['rt_time']]-spec2[method_dic[method]['rt_time']])<=0.1:
                node_a = f"{spec1[method_dic[method]['filename']]}_{spec1[method_dic[method]['scan']]}"
                node_b = f"{spec2[method_dic[method]['filename']]}_{spec2[method_dic[method]['scan']]}"
                G.add_edge(node_a,node_b)

    return G

def create_matching_network_across_files(cluster,method,tolerence):
    G = nx.Graph()

    for _, spec1 in cluster.iterrows():
        node_a =f"{spec1[method_dic[method]['filename']]}_{spec1[method_dic[method]['scan']]}"
        G.add_node(node_a, filename=spec1[method_dic[method]['filename']])
        for _, spec2 in cluster.iterrows():
            if (spec1[method_dic[method]['filename']] == spec2[method_dic[method]['filename']] and spec1[method_dic[method]['scan']] == spec2[method_dic[method]['scan']]):
                continue
            if abs(spec1['']):
                node_a =f"{spec1[method_dic[method]['filename']]}_{spec1[method_dic[method]['scan']]}"
                node_b =f"{spec2[method_dic[method]['filename']]}_{spec2[method_dic[method]['scan']]}"
                G.add_node(node_a,filename=spec1[method_dic[method]['filename']])
                G.add_node(node_b,filename=spec2[method_dic[method]['filename']])
                G.add_edge(node_a,node_b)
    return G


def calculate_max_component_per_file(G):

    # Find all connected components in the graph
    components = nx.connected_components(G)

    # Initialize a dictionary to hold the maximum component size for each file
    max_component_sizes = defaultdict(int)

    # Iterate through each component
    for component in components:
        # Create a temporary dictionary to count the number of nodes per file in this component
        file_counts = defaultdict(int)

        # Count nodes per file in the current component
        for node in component:
            filename = G.nodes[node]['filename']
            file_counts[filename] += 1

        # Update the max component size for each file encountered in this component
        for filename, count in file_counts.items():
            if count > max_component_sizes[filename]:
                max_component_sizes[filename] = count
    # Ensure that files represented by single nodes are accounted for
    for node in G.nodes:
        filename = G.nodes[node]['filename']
        if filename not in max_component_sizes:
            max_component_sizes[filename] = 1
        else:
            # Ensure there's at least a count of 1 for each file,
            # useful if a file's node(s) were not part of any component counted above
            max_component_sizes[filename] = max(max_component_sizes[filename], 1)

    return dict(max_component_sizes)

# Function to calculate purity for a cluster
def calculate_cluster_purity(cluster,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, method)

    max_component_sizes = calculate_max_component_per_file(G)
    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster[method_dic[method]['filename']])

    max_fraction = 0.0

    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        fraction = max_component_sizes[filename] /count
        max_fraction = max(fraction, max_fraction)

    return max_fraction

def calculate_cluster_purity_weighted_avg(task,method):
    _,cluster = task
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = optimized_create_matching_network(cluster, method)

    max_component_sizes = calculate_max_component_per_file(G)

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster[method_dic[method]['filename']])

    frequencies = []
    values = []
    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max_component_sizes[filename]
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    weighted_sum = sum(value * frequency for value, frequency in zip(values, frequencies))

    # Calculate the total frequency
    total_frequency = sum(frequencies)

    # Calculate the weighted average
    weighted_average = weighted_sum / total_frequency
    return weighted_average

def calculate_cluster_purity_avg(cluster,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster,method)

    # Get connected components
    max_component_sizes = calculate_max_component_per_file(G)
    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster[method_dic[method]['filename']])

    frequencies = []
    values = []
    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max_component_sizes[filename]
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    total_purity = sum(values)

    # Calculate the total frequency
    total_frequency = sum(frequencies)

    # Calculate the weighted average
    average = total_purity / total_frequency
    return average


def compare_scans(scan1, scan2, mass_tolerance=0.01, rt_tolerance=6000):
    mass_diff = abs(scan1[0] - scan2[0])
    rt_diff = abs(scan1[1] - scan2[1])
    return mass_diff <= mass_tolerance and rt_diff <= rt_tolerance

def assign_size_range_bin(size):
    if size == 1:
        return '1 to 1'  # Change here to ensure sorting works
    else:
        upper_limit = 2
        while size > upper_limit:
            upper_limit *= 2
        return f"{upper_limit//2+1} to {upper_limit}"


def apply_parallel_with_tqdm(groups, func, method):
    tasks = [(name, group) for name, group in groups]  # Ensure tasks are prepared as tuples
    partial_func = partial(func, method=method)

    with Pool(processes=60) as pool:
        result_list = []
        for result in tqdm(pool.imap(partial_func, tasks), total=len(tasks)):
            result_list.append(result)

    return result_list
def mscluster_purity(cluster_results):
    #handle the new version workflow filename issue
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzml')

    #return cluster_results.groupby('#ClusterIdx').progress_apply(lambda x: calculate_cluster_purity_weighted_avg(x, matching_pairs_set,'mscluster')), cluster_results.groupby('#ClusterIdx').size()
    return apply_parallel_with_tqdm(cluster_results.groupby('#ClusterIdx'),calculate_cluster_purity_weighted_avg,'mscluster'), cluster_results.groupby('#ClusterIdx').size()

def falcon_purity(cluster_results):
    #return cluster_results.groupby('cluster').progress_apply(lambda x: calculate_cluster_purity_weighted_avg(x,'falcon')), cluster_results.groupby('cluster').size()
    return apply_parallel_with_tqdm(cluster_results.groupby('cluster'), calculate_cluster_purity_weighted_avg, 'falcon'), cluster_results.groupby('cluster').size()
def maracluster_purity(cluster_results):
    #return cluster_results.groupby('cluster').progress_apply(lambda x: calculate_cluster_purity_weighted_avg(x, 'maracluster')), cluster_results.groupby('cluster').size()
    return apply_parallel_with_tqdm(cluster_results.groupby('cluster'), calculate_cluster_purity_weighted_avg,
                                    'maracluster'), cluster_results.groupby('cluster').size()
def generate_matching_pairs_set_file(folder_path,data_folder_path):
    matching_pairs_all_files = []
    for filename in tqdm(os.listdir(folder_path)):
        if filename.endswith('.mzML'):
            mzml_file = os.path.join(folder_path, filename)
            ms1_data = get_ms1_data(mzml_file)

            matching_pairs = []

            for scan1, scan2 in combinations(ms1_data, 2):
                if compare_scans(scan1, scan2):
                    matching_pairs.append((scan1[2], scan2[2]))

            matching_pairs_all_files.extend([(filename, pair[0], pair[1]) for pair in matching_pairs])

    matching_pairs_all_files = [(f'mzML/{filename}', scan_start, scan_end) for filename, scan_start, scan_end in
                                matching_pairs_all_files]
    matching_pairs_set = {(item[0], item[1], item[2]) for item in matching_pairs_all_files}
    with open(data_folder_path+'/matching_pairs_set.pkl', 'wb') as file:
        pickle.dump(matching_pairs_set, file)
    return matching_pairs_set

def mscluster_merge(cluster_results,correction = False):
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzml')
    if correction:
        cluster_results['#RetTime'] = cluster_results['#RetTime']


    # Create a merged dataset with left join to include all cluster results and match with database results

    return cluster_results
def falcon_merge(cluster_results):

    return cluster_results
def maracluster_merge(cluster_results,reference_data):
    # reference_data.drop('#ClusterIdx', axis=1, inplace=True)
    # reference_data['#Filename'] = reference_data['#Filename'].str.replace('input_spectra', 'mzml')
    # reference_data['#RetTime'] = reference_data['#RetTime']
    # cluster_results['filename'] = cluster_results['filename'].str.replace('data','mzml')
    # merged_data = pd.merge(cluster_results, reference_data, left_on=['filename', 'scan'], right_on=['#Filename', '#Scan'],how='inner')
    return cluster_results
def custom_sort(bin_range):
    if bin_range == '1 to 1':
        return 0, 1  # This ensures "1 to 1" has the smallest possible sort key
    else:
        # Extract the start and end numbers from the bin label and use them as sort keys
        start, end = map(int, bin_range.split(' to '))
        return start, end

def range_avg(cluster_purity, cluster_size):
    purity_by_cluster_size = pd.DataFrame({'ClusterSize': cluster_size, 'ClusterPurity': cluster_purity})
    purity_by_cluster_size['Bin'] = purity_by_cluster_size['ClusterSize'].apply(assign_size_range_bin)

    # Group by 'Bin' and calculate the mean purity for each bin
    average_purity_by_bin = purity_by_cluster_size.groupby('Bin')['ClusterPurity'].mean().reset_index()

    count_clusters_by_bin = purity_by_cluster_size.groupby('Bin')['ClusterSize'].count().reset_index()
    count_clusters_by_bin.rename(columns={'ClusterSize': 'ClusterCount'}, inplace=True)

    # Merge the two dataframes on 'Bin'
    merged_df = pd.merge(average_purity_by_bin, count_clusters_by_bin, on='Bin')

    # Sort by 'Bin'
    merged_df = merged_df.sort_values('Bin', key=lambda x: x.map(custom_sort))

    return merged_df

def calculate_weighted_average_purity(purity, cluster_size):
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

def purity_wapper(clustering_ressults,methods):
    if methods == 'mscluster':
        return mscluster_purity(clustering_ressults)
    elif methods == 'falcon':
        return falcon_purity(clustering_ressults)
    elif methods == 'maracluster':
        return maracluster_purity(clustering_ressults)

if __name__ == "__main__":
    tqdm.pandas()
    parser = argparse.ArgumentParser(description='Using ClassyFIre results to benchmark the network')
    parser.add_argument('-c', type=str, required=True, default="cluster_info.tsv", help='input clustering reuslts filename')
    parser.add_argument('-t',type=int, required=True,default = 109333, help='Number of MS/MS in the datasets')
    parser.add_argument('-methods', type=str, required=True, default="falcon", help='Clustering methods')
    args = parser.parse_args()
    clustering_filename = args.c
    methods = args.methods
    total_num_spec = args.t

    # Load cluster results
    clustering_file_path = '../data/'+methods+'/'+clustering_filename
    cluster_results = pd.read_csv(clustering_file_path, sep='\t')

    #calculate n10 number & average purity
    cluster_n10 = n10_wrapoer(cluster_results, total_num_spec,methods)
    print("N10 value for {} clustering results is {}".format(methods,cluster_n10))
    clustering_purity,clustering_size = purity_wapper(cluster_results,methods)
    avg_purity = sum(clustering_purity) / len(clustering_purity)
    print("Average Purity for {} is :{}".format(methods,avg_purity))