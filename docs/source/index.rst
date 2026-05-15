ePLACE Documentation
====================

**ePLACE** (environmental Phylogenetic Localisation and Clade Estimation) is a Python library for analyzing environmental DNA (eDNA) sequences through BLAST comparison and taxonomic classification.

.. image:: https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4
   :target: https://edwards.flinders.edu.au/
   :alt: Edwards Lab

.. image:: https://zenodo.org/badge/1120756704.svg
   :target: https://doi.org/10.5281/zenodo.18181123
   :alt: DOI

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

Overview
--------

ePLACE provides a comprehensive toolkit for environmental DNA analysis with the following capabilities:

* **NCBI Database Management**: Download and manage BLAST core_nt and MMseqs2 NT databases
* **FASTA File Processing**: Read and validate FASTA files
* **Sequence Search**: Run BLAST or MMseqs2 searches with configurable parameters
* **Result Filtering**: Filter BLAST results by identity and coverage thresholds
* **Taxonomic Analysis**: Extract representative sequences per taxonomic rank
* **Sequence Extraction**: Retrieve sequences from BLAST databases
* **Sequence Trimming**: Trim reference sequences to aligned regions based on BLAST coordinates
* **Multiple Sequence Alignment**: Align sequences using MAFFT with auto-orientation
* **Phylogenetic Trees**: Build and label phylogenetic trees using IQTree
* **Tree Relabeling**: Relabel existing trees with taxonomic names at different ranks
* **Results Summary Output**: Creates a tab separated output that summarises the per-sequence matches

Quick Start
-----------

Installation
~~~~~~~~~~~~

.. code-block:: bash
   # create and activate a mamba environment
   mamba create -yn eplace bioconda::blast bioconda::mmseqs2 bioconda::pytaxonkit bioconda::iqtree bioconda::mafft
   mamba activate eplace

   # Clone the repository
   git clone https://github.com/linsalrob/eplace.git
   cd eplace

   # Install the package
   pip install .

   # Or install in development mode
   pip install -e .[dev]

Basic Usage
~~~~~~~~~~~

.. code-block:: bash

   # Download NCBI BLAST core_nt database
   eplace download

   # Download MMseqs2 NT database
   eplace download --target mmseqs2

   # Run search analysis (BLAST by default)
   eplace search query.fasta output_dir

   # Run grouped analysis
   eplace grouped query.fasta output_dir --group-rank order
   
   # Relabel existing tree with taxonomic names
   eplace relabel blast_results.txt input.treefile output.treefile --rank genus

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   cli
   workflows

.. toctree::
   :maxdepth: 2
   :caption: Detailed Documentation

   ncbi_download
   blast_workflow
   alignment
   phylogenetic_trees

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api

.. toctree::
   :maxdepth: 1
   :caption: Additional Information

   contributing
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
