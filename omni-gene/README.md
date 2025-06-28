# Omni-Gene: An AI-Powered Multi-Omics Platform

Omni-Gene is a professional, command-line-driven bioinformatics platform designed for the end-to-end analysis of multi-omics datasets. It provides a fully automated workflow for data ingestion, quality control, multi-modal integration, and AI-powered discovery of novel biological insights.

The tool is designed for computational biologists and researchers who need a reproducible and robust way to analyze complex datasets combining genomics, transcriptomics, and proteomics.

## Core Features

- **Multi-Modal by Design:** Natively handles RNA, protein, and genomic variant data.
- **Packaged & Installable:** Built as a proper Python package with a simple command-line interface (`omnigene`).
- **Automated Workflow:** Executes a multi-step pipeline from raw data to final report with a single command.
-   **Intelligent Analysis:**
    -   Performs data-driven clustering to identify novel cell subtypes.
    -   Conducts differential analysis to find key features between conditions.
    -   Generates a final, publication-ready report archive.

---

## Installation

Due to the complexity of bioinformatics software, installation is best performed within a Linux environment. For Windows users, the **Windows Subsystem for Linux (WSL)** is the recommended and most reliable method.

### 1. (For Windows Users) Set Up WSL

If you are on Windows, you must complete this one-time setup first.

1.  **Install WSL:** Open **PowerShell as an Administrator** and run:
    ```powershell
    wsl --install
    ```
2.  **Restart your computer.**
3.  **Set up Ubuntu:** An Ubuntu terminal will open. Create a username and password when prompted.
4.  **Install Miniconda inside Ubuntu:** Open the Ubuntu terminal and run the following commands one by one:
    ```bash
    # Download the Miniconda installer
    wget [https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)

    # Run the installer. Accept all defaults by pressing ENTER and typing 'yes'.
    bash Miniconda3-latest-Linux-x86_64.sh

    # Close and re-open the Ubuntu terminal for changes to take effect.
    ```

### 2. Install Omni-Gene

These steps should be performed inside your **Ubuntu (WSL) terminal** or a native Linux/macOS terminal.

1.  **Clone the Project:**
    ```bash
    # Example for cloning to your desktop
    # Note: On WSL, your C: drive is at /mnt/c/
    cd /mnt/c/Users/YourUsername/Desktop/
    git clone <your-github-repo-url>
    cd omni-gene-project
    ```

2.  **Create and Activate Conda Environment:**
    This command reads the `environment.yaml` file and installs all dependencies in an isolated environment.
    ```bash
    conda env create -f environment.yaml
    conda activate omnigene_env
    ```

3.  **Install the Omni-Gene Package:**
    This command installs your tool in "editable" mode, so code changes are reflected immediately.
    ```bash
    pip install -e .
    ```

---

## How to Run

1.  **Prepare Your Data:**
    - Place your data files (e.g., `.h5ad`, `.csv`, `.vcf.gz`) in the appropriate subdirectories inside the `data/` folder.
    - Update your `metadata.csv` file with sample information.

2.  **Configure Your Analysis:**
    - Open and edit the main configuration file provided with the project to point to your data files and set analysis parameters.

3.  **Execute the Pipeline:**
    - Make sure your `omnigene_env` conda environment is active.
    - From the `omni-gene-project` directory, run the tool using the `omnigene` command:
    ```bash
    omnigene --config path/to/your/config.yaml
    ```

4.  **Check Your Results:**
    - A new timestamped folder will be created inside the `outputs/` directory. This folder will contain a `results.zip` file with all generated plots, tables, and reports.

