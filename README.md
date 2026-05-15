[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au/)
[![DOI](https://zenodo.org/badge/1120756704.svg)](https://doi.org/10.5281/zenodo.18181123)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/eplace)
[![Documentation Status](https://readthedocs.org/projects/eplace/badge/?version=latest)](https://eplace.readthedocs.io/en/latest/)

# eplace

ePLACE: environmental Phylogenetic Localisation and Clade Estimation

A Python library for analyzing environmental DNA (eDNA) sequences through BLAST comparison and taxonomic classification.

## Documentation

For all the features available, [please check out readthedocs](https://eplace.readthedocs.io/en/latest/)

## Features

- **NCBI Database Management**: Download and manage NCBI BLAST databases (core_nt)
- **FASTA File Processing**: Read and validate FASTA files
- **BLAST Search**: Run blastn searches with configurable parameters
- **Result Filtering**: Filter BLAST results by identity and coverage thresholds
- **Taxonomic Analysis**: Extract representative sequences per taxonomic rank
- **Sequence Extraction**: Retrieve sequences from BLAST databases
- **Sequence Trimming**: Trim reference sequences to aligned regions based on BLAST coordinates
- **Multiple Sequence Alignment**: Align sequences using MAFFT with auto-orientation
- **Phylogenetic Trees**: Build and label phylogenetic trees using IQTree
- **Tree Relabeling**: Relabel existing trees with taxonomic names at different ranks
- **Results Summary Output**: Creates a tab separated output that summarises the per-sequence matches.

## Installation

## conda: coming soon

## pip: the easy way

Using pip:

```bash
# create and activate a mamba environment
mamba create -yn eplace bioconda::blast bioconda::mmseqs2 bioconda::pytaxonkit bioconda::iqtree bioconda::mafft
mamba activate eplace

pip install eplace
```

### For development (not recommended)

```bash
# create and activate a mamba environment
mamba create -yn eplace bioconda::blast bioconda::mmseqs2 bioconda::pytaxonkit bioconda::iqtree bioconda::mafft
mamba activate eplace

# Clone the repository
git clone https://github.com/linsalrob/eplace.git
cd eplace

# Install the package
pip install -e .
```

After installation, the `eplace` command will be available in your environment.

## Requirements

- Python 3.8 or higher
- BLAST+ tools (blastn, blastdbcmd) must be installed separately:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install ncbi-blast+
  
  # macOS with Homebrew
  brew install blast
  ```
- TaxonKit (for taxonomy lookup):
  ```bash
  # Download from: https://github.com/shenwei356/taxonkit/releases
  # Or install via conda:
  conda install -c bioconda taxonkit
  ```
- MAFFT (optional, for sequence alignment):
  ```bash
  # Ubuntu/Debian
  sudo apt-get install mafft
  
  # macOS with Homebrew
  brew install mafft
  
  # Or via conda:
  conda install -c bioconda mafft
  ```
- IQTree (optional, for phylogenetic tree building):
  ```bash
  # Ubuntu/Debian
  sudo apt-get install iqtree
  
  # macOS with Homebrew
  brew install iqtree
  
  # Or via conda:
  conda install -c bioconda iqtree
  ```

## Quick Start

ePLACE provides a unified command-line interface with four main commands:

### 1. Download NCBI Database

```bash
# Download the BLAST core_nt database to default location
eplace download

# Download MMseqs2 NT database
eplace download --target mmseqs2

# Download both BLAST and MMseqs2 databases
eplace download --target both

# Download MMseqs2 NT and add taxonomy sidecar files
eplace download --target mmseqs2 --add-taxonomy --ncbi-taxonomy /path/to/ncbi/taxonomy/current

# Force redownload even if database exists
eplace download --force
```

### 2. Run Individual Search Analysis

Run sequence search and build one phylogenetic tree per query sequence:

```bash
# Basic usage with default parameters
eplace search query.fasta output_dir

# With custom parameters
eplace search query.fasta output_dir \
    --rank genus \
    --min-identity 95 \
    --min-coverage 85 \
    --num-threads 4

# Skip alignment and tree building (search and extraction only)
eplace search query.fasta output_dir --skip-alignment

# Show help
eplace search --help
```

### 3. Run Grouped Search Analysis

Run sequence search and group queries by taxonomic rank for joint phylogenetic analysis:

```bash
# Basic usage (group by class, default)
eplace grouped query.fasta output_dir

# Group by different taxonomic rank
eplace grouped query.fasta output_dir --group-rank order

# Specify both representative and grouping ranks
eplace grouped query.fasta output_dir --rank genus --group-rank family

# Show help
eplace grouped --help
```

### 4. Relabel Phylogenetic Trees

Relabel an existing phylogenetic tree with taxonomic names from BLAST results:

```bash
# Relabel tree with genus names
eplace relabel blast_results.txt input.treefile output.treefile --rank genus

# Relabel tree with species names (genus + species)
eplace relabel blast_results.txt input.treefile output.treefile --rank species

# Relabel tree with family names
eplace relabel blast_results.txt input.treefile output.treefile --rank family

# Show help
eplace relabel --help
```

### Using the Library API

You can also use ePLACE as a Python library:

```python
from eplace_lib import setup_ncbi_database

# Download the core_nt database
success, message = setup_ncbi_database()
```

```python
from pathlib import Path
from eplace_lib import run_blast_search, process_blast_results_for_taxonomy

# Run BLAST search with filtering
success, filtered_hits = run_blast_search(
    query_fasta=Path("query.fasta"),
    output_file=Path("blast_results.txt"),
    min_identity=90.0,    # 90% identity threshold
    min_coverage=80.0     # 80% query coverage threshold
)

# Extract representative sequences by taxonomic rank
results = process_blast_results_for_taxonomy(
    blast_hits=filtered_hits,
    output_dir=Path("output"),
    rank="species"  # Options: phylum, class, order, family, genus, species
)
```

## Command-Line Interface

The `eplace` command provides four subcommands:

### eplace download

Download and setup the NCBI core_nt BLAST database and/or MMseqs2 NT database.

**Usage:**
```bash
eplace download [--target {blast,mmseqs2,both}] [--force] [MMSEQS_OPTIONS]
```

**Options:**
- `--target {blast,mmseqs2,both}`: Which backend database(s) to download (default: `blast`)
- `--force`: Force redownload even if database exists
- `--mmseqs-db-dir PATH`: MMseqs2 database root directory (default: `$MMSEQS_DB_DIR`, then `$MMSEQS2DB`, or `~/mmseqs2db`)
- `--mmseqs-threads INT`: Threads for MMseqs2 download/taxonomy steps
- `--add-taxonomy`: Build MMseqs2 taxonomy sidecar files (`mmseqs createtaxdb`) after NT download
- `--ncbi-taxonomy PATH`: Path to NCBI taxonomy dump containing `nodes.dmp`, `names.dmp`, `merged.dmp` (required with `--add-taxonomy`)
- `--acc2taxid-dir PATH`: Path to accession2taxid files (defaults to `$ACC2TAXID_DIR` or `<ncbi-taxonomy>/accession2taxid`)
- `--taxonomy-workdir PATH`: Working directory for MMseqs taxonomy mapping files
- `--skip-memory-check`: Skip MMseqs2 RAM preflight checks

**Notes:**
- Database will be stored in `$BLASTDB` if set, otherwise `~/blastdb`
- MMseqs2 database is stored in `$MMSEQS_DB_DIR`, then `$MMSEQS2DB`, or `~/mmseqs2db`
- MMseqs2 NT download typically requires at least 64 GiB RAM
- MMseqs2 taxonomy integration typically requires at least 128 GiB RAM
- Downloads are large and may take time

### eplace search

Run sequence search with individual taxonomy analysis. Creates one phylogenetic tree per query sequence.

**Usage:**
```bash
eplace search QUERY_FASTA OUTPUT_DIR [OPTIONS]
```

**Required Arguments:**
- `QUERY_FASTA`: Path to query FASTA file
- `OUTPUT_DIR`: Output directory for results

**Optional Arguments:**
- `--rank {phylum,class,order,family,genus,species}`: Taxonomic rank for representative selection (default: genus)
- `--tree-label-rank {phylum,class,order,family,genus,species}`: Taxonomic rank for tree labeling (default: genus)
- `--min-identity FLOAT`: Minimum percent identity for BLAST hits (default: 90.0)
- `--min-coverage FLOAT`: Minimum query coverage percentage (default: 80.0)
- `--database NAME`: BLAST database name (default: core_nt)
- `--blastdb-path PATH`: Path to BLAST database directory
- `--num-threads INT`: Number of threads for BLAST and alignment (default: 1)
- `--search-tool {blast,mmseqs2}`: Sequence search backend to use (default: blast)
- `--blast-db-source LABEL`: BLAST database provenance label for `search_metadata.json` (defaults to `--database`)
- `--mmseqs-database NAME`: MMseqs2 database name (defaults to `--database` when using MMseqs2)
- `--mmseqs-db-path PATH`: MMseqs2 database directory (defaults to `$MMSEQS_DB_DIR`, then `$MMSEQS2DB`, then `~/mmseqs2db`)
- `--mmseqs-db-source LABEL`: MMseqs2 provenance label recorded in `search_metadata.json`
- `--mmseqs-sensitivity FLOAT`: MMseqs2 sensitivity, range 1-7.5 (default: 5.7)
- `--mmseqs-search-type INT`: MMseqs2 search type for `easy-search` (default: 3)
- `--overwrite-existing-blast`: Overwrite existing BLAST results
- `--skip-alignment`: Skip alignment and tree building steps
- `--output-classification PATH`: Path to output classification TSV file

### eplace grouped

Run BLAST search with grouped taxonomy analysis. Groups queries by taxonomic rank and creates one phylogenetic tree per group.

**Usage:**
```bash
eplace grouped QUERY_FASTA OUTPUT_DIR [OPTIONS]
```

**Required Arguments:**
- `QUERY_FASTA`: Path to query FASTA file
- `OUTPUT_DIR`: Output directory for results

**Optional Arguments:**
- `--rank {phylum,class,order,family,genus,species}`: Taxonomic rank for representative selection (default: genus)
- `--group-rank {phylum,class,order,family,genus,species}`: Taxonomic rank for grouping sequences (default: class)
- `--tree-label-rank {phylum,class,order,family,genus,species}`: Taxonomic rank for tree labeling (default: genus)
- `--combined-tree-label-rank {phylum,class,order,family,genus,species}`: Taxonomic rank for labeling the combined tree. If not provided, the combined tree will not be built (optional)
- `--min-identity FLOAT`: Minimum percent identity for BLAST hits (default: 90.0)
- `--min-coverage FLOAT`: Minimum query coverage percentage (default: 80.0)
- `--database NAME`: BLAST database name (default: core_nt)
- `--blastdb-path PATH`: Path to BLAST database directory
- `--num-threads INT`: Number of threads for BLAST and alignment (default: 1)
- `--search-tool {blast,mmseqs2}`: Sequence search backend to use (default: blast)
- `--blast-db-source LABEL`: BLAST database provenance label for `search_metadata.json` (defaults to `--database`)
- `--mmseqs-database NAME`: MMseqs2 database name (defaults to `--database` when using MMseqs2)
- `--mmseqs-db-path PATH`: MMseqs2 database directory (defaults to `$MMSEQS_DB_DIR`, then `$MMSEQS2DB`, then `~/mmseqs2db`)
- `--mmseqs-db-source LABEL`: MMseqs2 provenance label recorded in `search_metadata.json`
- `--mmseqs-sensitivity FLOAT`: MMseqs2 sensitivity, range 1-7.5 (default: 5.7)
- `--mmseqs-search-type INT`: MMseqs2 search type for `easy-search` (default: 3)
- `--overwrite-existing-blast`: Overwrite existing BLAST results
- `--skip-alignment`: Skip alignment and tree building steps
- `--alignment-tolerance INT`: Maximum coordinate difference for alignment consistency (default: 50)
- `--output-classification PATH`: Path to output classification TSV file

**Note:** The grouped workflow creates individual trees for each taxonomic group. Optionally, you can also create a combined tree from all groups by specifying `--combined-tree-label-rank`. The combined tree includes representatives from all taxonomic groups and can be very time-consuming to build with large datasets, so it is only built when explicitly requested.

### eplace relabel

Relabel a phylogenetic tree with taxonomic names from BLAST results. This is useful when you have an existing tree and want to replace sequence IDs with taxonomic names, or when you want to relabel a tree at a different taxonomic rank.

**Usage:**
```bash
eplace relabel BLAST_OUTPUT TREE_FILE OUTPUT_TREE [OPTIONS]
```

**Required Arguments:**
- `BLAST_OUTPUT`: Path to BLAST output file (tabular format with taxonomy)
- `TREE_FILE`: Path to input tree file (Newick format)
- `OUTPUT_TREE`: Path to output relabeled tree file

**Optional Arguments:**
- `--rank {phylum,class,order,family,genus,species}`: Taxonomic rank for tree labeling (default: genus)
- `--blastdb-path PATH`: Path to BLAST database directory (optional, not required for relabeling)

**Key Features:**
- Supports all standard taxonomic ranks from phylum to species
- Handles species names as "genus species" format for binomial nomenclature
- Preserves tree topology while updating labels
- Works with Newick format trees
- Handles reversed sequences (with _R_ prefix from MAFFT)
- Cleans labels for Newick format compatibility

**Examples:**
```bash
# Relabel tree with genus names (default)
eplace relabel blast_results.txt input.treefile output_labeled.treefile

# Relabel tree with species names (genus + species binomial)
eplace relabel blast_results.txt input.treefile output_species.treefile --rank species

# Relabel tree with family names
eplace relabel blast_results.txt input.treefile output_family.treefile --rank family

# Relabel tree using custom BLAST database location
eplace relabel blast_results.txt input.treefile output.treefile --rank genus --blastdb-path /path/to/blastdb
```

**Use Cases:**
- Re-label trees at different taxonomic ranks without rebuilding
- Add taxonomic labels to trees from external phylogenetic tools
- Create multiple versions of the same tree with different label granularity
- Update tree labels when taxonomy information changes

## Documentation

Full documentation is available at [Read the Docs](https://eplace.readthedocs.io/).

- [Installation Guide](https://eplace.readthedocs.io/en/latest/installation.html) - Complete installation instructions
- [Quick Start Guide](https://eplace.readthedocs.io/en/latest/quickstart.html) - Get started quickly
- [Command-Line Interface](https://eplace.readthedocs.io/en/latest/cli.html) - Complete CLI reference
- [Workflows](https://eplace.readthedocs.io/en/latest/workflows.html) - Individual and grouped workflow details
- [API Reference](https://eplace.readthedocs.io/en/latest/api.html) - Python API documentation
- [NCBI Database Download](https://eplace.readthedocs.io/en/latest/ncbi_download.html) - Database management guide
- [BLAST Workflow](https://eplace.readthedocs.io/en/latest/blast_workflow.html) - BLAST analysis guide
- [Alignment](https://eplace.readthedocs.io/en/latest/alignment.html) - Sequence alignment documentation
- [Phylogenetic Trees](https://eplace.readthedocs.io/en/latest/phylogenetic_trees.html) - Tree building guide
- [Contributing](https://eplace.readthedocs.io/en/latest/contributing.html) - Contribution guidelines

### Local Documentation

You can also build the documentation locally:

```bash
cd docs
make html
# Open docs/build/html/index.html in your browser
```

## Workflow Comparison

### Individual Workflow (`eplace search`)

The individual workflow processes each query sequence independently:
- Creates one output directory per query sequence
- Extracts representative sequences for each query at the specified taxonomic rank
- Builds one multiple sequence alignment per query
- Creates one phylogenetic tree per query

**Use when:** You want to analyze each query sequence in its own phylogenetic context.

### Grouped Workflow (`eplace grouped`)

The grouped workflow combines queries by taxonomic classification:
- Groups all queries that match to the same taxonomic rank (e.g., class, order)
- Creates one FASTA file per group containing all queries and unique reference sequences
- Removes redundant reference sequences within each group
- Builds one alignment and phylogenetic tree per group (instead of per query)

**Use when:** You want to analyze multiple related queries together in a single phylogenetic context.

### Examples

```bash
# Group queries by class (default)
eplace grouped query.fasta output_dir

# Group by a different taxonomic rank
eplace grouped query.fasta output_dir --group-rank order

# Specify both representative rank and grouping rank
eplace grouped query.fasta output_dir --rank genus --group-rank family
```

## Testing

Run the test suite:

```bash
# Run all tests
pytest tests/ -v

# Run specific test modules
pytest tests/test_blast_analysis.py -v
pytest tests/test_taxonomy.py -v
pytest tests/test_workflow.py -v

# Run with coverage
pytest tests/ --cov=eplace_lib --cov-report=html
```

## Project Structure

```
eplace/
├── src/
│   └── eplace_lib/
│       ├── __init__.py
│       ├── blast_analysis.py    # BLAST operations
│       ├── ncbi_download.py     # Database management
│       ├── sequences.py         # Sequence analysis utilities
│       └── taxonomy.py          # Taxonomy extraction
├── tests/
│   ├── test_blast_analysis.py
│   ├── test_ncbi_download.py
│   ├── test_taxonomy.py
│   └── test_workflow.py
├── examples/
│   ├── blast_workflow_example.py
│   └── download_ncbi_example.py
├── docs/
│   ├── blast_workflow.md
│   └── ncbi_download.md
└── pyproject.toml
```

## Workflow Overview

1. **Download Database**: Use `setup_ncbi_database()` to download NCBI core_nt database
2. **Prepare Query**: Create a FASTA file with your query sequences
3. **Run BLAST**: Use `run_blast_search()` to search against the database
4. **Filter Results**: Automatically filter by identity and coverage thresholds
5. **Extract Representatives**: Select representative sequences per taxonomic rank
6. **Trim Sequences**: Extract aligned regions from reference sequences based on BLAST coordinates
7. **Align Sequences**: Use MAFFT to align query with trimmed reference sequences (optional)
8. **Build Tree**: Build phylogenetic tree using IQTree with taxonomic labels (optional)
9. **Output**: Get FASTA files, alignments, and trees (one set per query)

### Grouped Workflow Overview

The grouped workflow adds an additional step:
1-5. Same as standard workflow through representative extraction
6. **Group by Rank**: Group all queries by specified taxonomic rank (e.g., class)
7. **Create Grouped FASTA**: Combine all queries and unique references for each group
8. **Trim Sequences**: Trim references to aligned regions
9. **Check Consistency**: Verify BLAST hits align to similar locations on references
10. **Align and Build Trees**: Create one alignment and tree per taxonomic group
11. **Build Combined Tree**: Create a combined tree from all groups with all queries and representatives

## Output Structure

### Standard Workflow Output

```
output_dir/
├── blast_results.txt              # Raw BLAST results
├── blast_results_annotated.txt    # BLAST results with taxonomic annotations
├── query1_id/
│   ├── query1_id_representatives.fasta          # Representative sequences
│   ├── query1_id_with_query.fasta              # Query + representatives
│   ├── query1_id_trimmed.fasta                 # Trimmed to aligned regions
│   ├── query1_id_aligned.fasta                 # Multiple sequence alignment
│   ├── query1_id_tree.treefile                 # Phylogenetic tree
│   ├── query1_id_tree_labeled.treefile         # Tree with taxonomic labels
│   └── query1_id_tree.* (other IQTree files)
├── query2_id/
│   └── query2_id_representatives.fasta
└── ...
```

### Grouped Workflow Output

```
output_dir/
├── blast_results.txt              # Raw BLAST results
├── blast_results_annotated.txt    # BLAST results with taxonomic annotations
├── query1_id/                     # Per-query representative sequences (from step 5)
│   └── query1_id_representatives.fasta
├── query2_id/
│   └── query2_id_representatives.fasta
├── Taxonomic_Group_1/             # One directory per taxonomic group
│   ├── Taxonomic_Group_1_combined.fasta        # All queries + unique references
│   ├── Taxonomic_Group_1_trimmed.fasta         # Trimmed to aligned regions
│   ├── Taxonomic_Group_1_aligned.fasta         # Multiple sequence alignment
│   ├── Taxonomic_Group_1_tree.treefile         # Phylogenetic tree
│   ├── Taxonomic_Group_1_tree_labeled.treefile # Tree with taxonomic labels
│   └── Taxonomic_Group_1_tree.* (other IQTree files)
├── Taxonomic_Group_2/
│   └── ...
├── combined_all_groups_trimmed.fasta           # Combined alignment of all groups
├── combined_all_groups_aligned.fasta           # Multiple sequence alignment
├── combined_all_groups_tree.treefile           # Combined phylogenetic tree
├── combined_all_groups_tree_labeled.treefile   # Combined tree with taxonomic labels
└── ...
```

## License

MIT License - See LICENSE file for details

## Authors

- Rob Edwards (raedwards@gmail.com)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use ePLACE in your research, please cite:

```
Edwards, R. (2024). ePLACE: environmental Phylogenetic Localisation and Clade Estimation.
GitHub repository: https://github.com/linsalrob/eplace
```

## Support

For issues, questions, or suggestions, please open an issue on GitHub:
https://github.com/linsalrob/eplace/issues
