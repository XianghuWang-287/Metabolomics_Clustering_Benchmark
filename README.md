# Metabolomics_Clustering_Benchmark
Git Hub repo for the paper "Evaluation of Tandem Mass Spectra Clustering Tools for Metabolomics". Including tools to benchmark clusterng tools for metabolomics datasets

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
source Clustering_benchmarking/bin/activate
pip install -r requirements.txt
```

## Dataset Preparation 📊

All datasets used in this project are publicly available and can be downloaded from [ProteomeXchange](http://www.proteomexchange.org/), a globally coordinated initiative to facilitate the exchange and dissemination of proteomics data.

### Downloading Datasets from ProteomeXchange

1. **Visit the ProteomeXchange website:**
   - Navigate to [www.proteomexchange.org](http://www.proteomexchange.org/).

2. **Search for Datasets:**
   - Use the search bar to find datasets relevant to metabolomics clustering.
   - You can search by keywords, dataset identifiers, or browse through the available datasets.

3. **Download Raw Data Files:**
   - Select the datasets you wish to use.
   - Download the raw data files, which are typically in vendor-specific formats (e.g., Thermo `.raw` files).

### Converting Raw Data to mzML or MGF Format

To process the raw data files, they need to be converted into standard formats like **mzML** or **MGF**. We recommend using either **msconvert** or **ThermoRawFileParser** for this purpose.

---

#### **Option 1: Using msconvert**

**msconvert** is a tool from the ProteoWizard toolkit that converts vendor-specific mass spectrometry data formats into open formats.

- **Website:** [ProteoWizard Downloads](http://proteowizard.sourceforge.net/download.html)

##### Installation

1. **Download ProteoWizard:**
   - Visit the [ProteoWizard download page](http://proteowizard.sourceforge.net/download.html).
   - Choose the appropriate installer for your operating system (Windows or macOS).

2. **Install the Software:**
   - Run the installer and follow the on-screen instructions.

##### Conversion Instructions

1. **Open a Terminal or Command Prompt.**

2. **Navigate to Your Data Directory:**

   Use `cd` to change to the directory containing your raw data files.

3. **Convert to mzML Format:**

   ```bash
   msconvert your_file.raw --mzML
    ```
#### **Option 2: Using ThermoRawFileParser**

**ThermoRawFileParser** is a cross-platform tool for converting Thermo RAW files to mzML or MGF formats without the need for vendor software.

- **Website:** [ThermoRawFileParser GitHub Repository](https://github.com/compomics/ThermoRawFileParser)

##### Installation

1. **Download the Latest Release:**

   - Go to the [ThermoRawFileParser releases page](https://github.com/compomics/ThermoRawFileParser/releases).
   - Download the appropriate version for your operating system.

2. **Extract the Archive:**

   - Unzip the downloaded file to a directory of your choice.

##### Conversion Instructions

1. **Open a Terminal or Command Prompt.**

2. **Navigate to the ThermoRawFileParser Directory:**

   ```bash
   cd /path/to/ThermoRawFileParser
   ```
##### Conversion Instructions

3. **Convert to mzML Format:**
   - On Windows:

   ```bash
   ThermoRawFileParser.exe -i "C:\path\to\your_file.raw" -o "C:\path\to\output\folder" -f=1
   ```
   
   - On macOS/Linux (using Mono):
   
   ```bash
   mono ThermoRawFileParser.exe -i /path/to/your_file.raw -o /path/to/output/folder -f=1
   ```
   
4. **Convert to MGF Format:**
   - On Windows:
   
   ```bash
   ThermoRawFileParser.exe -i "C:\path\to\your_file.raw" -o "C:\path\to\output\folder" -f=0
   ```
   
   - On macOS/Linux (using Mono):
   
   ```bash
   mono ThermoRawFileParser.exe -i /path/to/your_file.raw -o /path/to/output/folder -f=0
   ```
   
5. **Additional Options:**
   -For more advanced usage, refer to the [ThermoRawFileParser documentation](https://github.com/compomics/ThermoRawFileParser).

## Post-processing of Clustering Tool Results 🛠️

After running the clustering tools, you'll need to prepare the results for analysis. This section provides step-by-step instructions for post-processing the results from **msCluster**, **Falcon**, and **MaRaCluster** to ensure they have consistent column names and formats required for downstream analysis.

---

### msCluster

For **msCluster**, you can use any processing protocol as long as the output file has the following column names:

- `'filename'`: `#Filename`
- `'scan'`: `#Scan`
- `'mass'`: `#ParentMass`
- `'rt_time'`: `#RetTime`
- `'cluster'`: `#ClusterIdx`

An example of the processed data can be found in the demo data folder:

- `./data/mscluster/mscluster_cluster_info.tsv`


#### Recommended Processing Method

The easiest way to process msCluster results is by using the **GNPS 2.0 Classical Networking Workflow**. Follow the steps below:

1. **Access the GNPS 2.0 Workflow:**

   - Navigate to the [GNPS 2.0 Classical Networking Workflow](https://gnps2.org/workflowinput?workflowname=classical_networking_workflow).

2. **Upload Your Data:**

   - Submit your original data files (e.g., mzML files) to the workflow.

3. **Run the Workflow:**

   - Configure the workflow parameters as needed and start the job.

4. **Retrieve the Processed Results:**

   - After the workflow completes, download the processed results from:

     ```
     nf_output/clustering/clusterinfo.tsv
     ```

   - This file will have the required column names and is ready for downstream analysis.

---

### Falcon

To process the results from **Falcon**, follow these steps:

1. **Run Falcon Clustering:**

   - Follow the official instructions provided on the [Falcon GitHub page](https://github.com/bittremieux/falcon) to perform clustering.

2. **Prepare the Results for Analysis:**

   - Use the provided script in this project to summarize and format the Falcon results.

   - **Script Location:**

     ```
     src/utility/summarize_results
     ```

   - **Command to Run the Script:**

     ```bash
     python3 src/utility/summarize_results ./falcon.csv output_summary
     ```

     - Replace `./falcon.csv` with the path to your Falcon results file.
     - `output_summary` is the desired name for the processed output file.

   - **Note:**

     - The script will format the results and ensure they have the necessary columns for analysis.

---

### MaRaCluster

Processing **MaRaCluster** results requires merging retention time information with the clustering results. Follow these steps:

1. **Run MaRaCluster Clustering:**

   - Follow the official instructions on the [MaRaCluster GitHub page](https://github.com/statisticalbiotechnology/maracluster) to perform clustering.

2. **Retrieve Retention Time Information:**

   - Use the **GNPS 2.0 PerScanSummarizer Workflow** to extract retention time information.

   - **Access the Workflow:**

     - Navigate to the [GNPS 2.0 PerScanSummarizer](https://gnps2.org/workflowinput?workflowname=PerScanSummarizer).

   - **Upload Your Data:**

     - Submit your original data files (e.g., mzML files) to the workflow.

   - **Run the Workflow:**

     - Configure the workflow parameters as needed and start the job.

   - **Download the Retention Time Results:**

     - After the workflow completes, download the results file:

       ```
       nf_output/results.csv
       ```

3. **Merge Retention Time with Clustering Results:**

   - Use the provided script in this project to merge the retention time information with the MaRaCluster results.

   - **Script Location:**

     ```
     src/utility/maracluster_processing
     ```

   - **Modify the Script Inputs:**

     - Open the script and set the following variables:

       - **MaRaCluster Results Path:**

         ```python
         input_file = 'path/to/maracluster_results.tsv'
         ```

       - **Retention Time Workflow Results Path:**

         ```python
         input_workflow_reuslt = 'path/to/nf_output/results.csv'
         ```

       - **Output File Path:**

         ```python
         output_file = 'path/to/output_file.tsv'
         ```

   - **Run the Script:**

     ```bash
     python3 src/utility/maracluster_processing
     ```

   - **Ensure Required Columns:**

     - The final output file should have at least the following columns:

       - `'filename'`: `filename`
       - `'scan'`: `scan`
       - `'mass'`: `precursor_mz`
       - `'rt_time'`: `retention_time`
       - `'cluster'`: `cluster`

---

## Running the Benchmarking Script 📈

Once you have post-processed the clustering results from **msCluster**, **Falcon**, and **MaRaCluster**, you can benchmark their performance using the provided benchmarking script. This section outlines how to execute the benchmarking script and explains the required parameters.

### Script Location

The benchmarking script is located in the `src` directory of the project:

```
src/Clustering_benchmark_MS_RT.py
```

### Prerequisites

Before running the benchmarking script, ensure that you have completed the following:

- **Post-processing Clustering Results:** Ensure that the clustering results from **msCluster**, **Falcon**, and **MaRaCluster** have been post-processed and are in the required format as described in the [Post-processing of Clustering Tool Results](#post-processing-of-clustering-tool-results) section.
- **Python Environment:** Ensure that your virtual environment is activated and all necessary dependencies are installed (as per the [Installation Guide](#installation-guide-🛠️)).

### Usage

To run the benchmarking script, use the following command structure:

```bash
python3 src/Clustering_benchmark_MS_RT.py -c <cluster_info_file> -t <number_of_msms> -methods <clustering_method> [-tol <tolerance>]
```

### Parameters

The benchmarking script accepts several command-line arguments to customize its behavior. Below is a detailed explanation of each parameter:

| Parameter | Type   | Required | Default    | Description                                           |
|-----------|--------|----------|------------|-------------------------------------------------------|
| `-c`      | `str`  | Yes      | `cluster_info.tsv` | **Input Clustering Results Filename:** Path to the clustering results file (e.g., `./data/mscluster/mscluster_cluster_info.tsv`). |
| `-t`      | `int`  | Yes      | `109333`   | **Number of MS/MS in the Datasets:** Total number of MS/MS spectra in your dataset. |
| `-methods`| `str`  | Yes      | `falcon`   | **Clustering Methods:** Specifies the clustering method to benchmark (e.g., `falcon`, `mscluster`, `maracluster`). |
| `-tol`    | `float`| No       | `0.1`      | **Tolerance for the MS-RT Window:** Tolerance value for the mass-to-retention time window used in benchmarking. |

### Example Commands

Below are example commands demonstrating how to run the benchmarking script for different clustering methods.

#### Benchmarking Falcon Results

```bash
python3 src/Clustering_benchmark_MS_RT.py -c ./data/falcon/falcon_cluster_info.tsv -t 109333 -methods falcon -tol 0.1
```

#### Benchmarking msCluster Results

```bash
python3 src/Clustering_benchmark_MS_RT.py -c ./data/mscluster/mscluster_cluster_info.tsv -t 109333 -methods mscluster -tol 0.1
```

#### Benchmarking MaRaCluster Results

```bash
python3 src/Clustering_benchmark_MS_RT.py -c ./data/maracluster/maracluster_cluster_info.tsv -t 109333 -methods maracluster -tol 0.1
```

### Parameter Descriptions

- **Input Clustering Results Filename (`-c`):**
  - **Description:** Path to the post-processed clustering results file.
  - **Example:** `./data/falcon/falcon_cluster_info.tsv`

- **Number of MS/MS in the Datasets (`-t`):**
  - **Description:** Total number of MS/MS spectra in your dataset. This value is used for normalization and benchmarking purposes.
  - **Example:** `109333`

- **Clustering Methods (`-methods`):**
  - **Description:** Specifies which clustering method's results you want to benchmark. Supported methods include `falcon`, `mscluster`, and `maracluster`.
  - **Example:** `falcon`

- **Tolerance for the MS-RT Window (`-tol`):**
  - **Description:** Defines the tolerance for the mass-to-retention time window when benchmarking. This parameter adjusts the strictness of matching between clustered results and reference data.
  - **Default Value:** `0.1`
  - **Example:** `0.05`

### Notes

- **Consistency Across Methods:** Ensure that the `-c` parameter points to the correctly post-processed clustering results corresponding to the specified `-methods` parameter.
- **Tolerance Adjustment:** You may need to adjust the `-tol` parameter based on the characteristics of your dataset to achieve optimal benchmarking results.
- **Multiple Methods:** To benchmark multiple clustering methods, run the script separately for each method with the corresponding `-methods` parameter.

### Troubleshooting

- **Missing Dependencies:** If you encounter import errors, ensure that all required Python packages are installed in your virtual environment.
- **Incorrect File Paths:** Verify that the paths provided to the `-c` parameter are correct and that the files exist.
- **Insufficient Permissions:** Ensure you have the necessary permissions to read input files and write output files in the specified directories.

---






