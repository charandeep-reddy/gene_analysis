# Gene Language Analytics Web Application

This is a web application for analyzing gene sequences using various computational methods.

## Features

- Sequence Alignment Analysis
- k-mer Frequency Analysis
- Motif and Pattern Detection
- Support for both text input and file upload
- Beautiful and responsive user interface

## Setup Instructions

1. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install the required packages:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
python app.py
```

4. Open your web browser and navigate to `http://localhost:5000`

## Usage

1. Enter or upload your gene sequences:
   - You can either paste the sequences directly into the text areas
   - Or upload text/FASTA files containing the sequences

2. Select an analysis method:
   - Sequence Alignment Analysis: Compares two sequences and shows their alignment
   - k-mer Frequency Analysis: Analyzes the frequency of k-length subsequences
   - Motif Detection: Identifies recurring patterns in the sequences

3. Click "Analyze" to process the sequences and view the results

## File Format Support

- Plain text files (.txt)
- FASTA format files (.fasta, .fna)
- Direct text input

## FASTA Files Included

The application comes with several FASTA (.fna) files that can be used for analysis:

- BRACA1.fna - Breast Cancer 1 gene sequence
- AOPE.fna - Apolipoprotein E gene sequence
- TP53.fna - Tumor protein p53 gene sequence
- CFTR.fna - Cystic Fibrosis Transmembrane Conductance Regulator gene sequence
- EGFR.fna - Epidermal Growth Factor Receptor gene sequence

These files are located in the root directory and can be used directly with the file upload feature.

## Requirements

See `requirements.txt` for a complete list of dependencies.