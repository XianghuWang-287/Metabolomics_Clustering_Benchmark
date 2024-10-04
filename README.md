# Metabolomics_Clustering_Benchmark
Git Hub repo for the paper "Evaluation of Tandem Mass Spectra Clustering Tools for Metabolomics". Including tools to benchmark clusteirng tools for metabolomics datasets

## Project Structure 📂

Below is the hierarchical outline of the key directories and files within this project:

```plaintext
📦 project_directory
 ┣ 📂 data
 ┃ ┣ 📂 mscluster
 ┃ ┣ 📂 falcon
 ┃ ┣ 📂 maracluster
 ┃ ┗ 📂 DB-Search
 ┣ 📂 src
 ┃ ┣ 📂 cluster_benchmark_MS_RT
 ┃ ┣ 📂 clustering_eval_compare_DB
 ┃ ┗ 📂 Plot
 ┗ 📂 results
   ┣ 📂 MS-RT
   ┗ 📂 DB-Search
📜 requirement.txt
📜 README.md
```

## Installation Guide 🛠️

Before running the project, ensure you have Python installed on your system. This project requires Python 3.8 or newer. You can download Python from [python.org](https://www.python.org/downloads/).

### Setting Up a Virtual Environment

It's recommended to use a virtual environment for Python projects. This keeps dependencies required by different projects separate by creating isolated environments for them. Here's how you can set it up:

1. Open a terminal or command prompt.
2. Navigate to the project's root directory.
3. Run the following command to create a virtual environment named `Clustering_benchmarking`:

```bash
python3 -m venv Clustering_benchmarking
source MN_benchmarking/bin/activate
pip install -r requirements.txt
```

## Prepare the Testing Dataset 📚

Before you can start the benchmarking process, you need to prepare the dataset that will be used for testing. This involves either downloading a pre-existing dataset or generating your own raw data in the required format. Follow the steps below to prepare your testing dataset:

### Downloading the Pre-existing Dataset