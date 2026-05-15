Quick Start Guide
=================

This guide will help you get started with ePLACE quickly.

Prerequisites
-------------

Before you begin, ensure you have:

1. Installed ePLACE and its dependencies (see :doc:`installation`)
2. BLAST+ tools installed (``blastn``, ``blastdbcmd``)
3. MMseqs2 installed (if you plan to use ``--search-tool mmseqs2``)
4. TaxonKit installed
5. (Optional) MAFFT and IQTree for alignment and phylogenetic analysis

Your First Analysis
-------------------

Step 1: Download Search Databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, download the database backend(s) you want to use:

.. code-block:: bash

   # Download BLAST core_nt database (default)
   eplace download

   # Download MMseqs2 NT database
   eplace download --target mmseqs2

   # Download both BLAST and MMseqs2 databases
   eplace download --target both

By default, BLAST core_nt is downloaded to your BLASTDB location
(``$BLASTDB`` or ``~/blastdb``). MMseqs2 databases use ``$MMSEQS_DB_DIR``,
then ``$MMSEQS2DB``, then ``~/mmseqs2db``.

.. note::
   These are large downloads and may take some time. You only need to do this once per backend.

Step 2: Prepare Your Query Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a FASTA file with your query sequences:

.. code-block:: text

   >query_sequence_1
   ATGCATGCATGCATGCATGCATGC
   >query_sequence_2
   GCTAGCTAGCTAGCTAGCTAGCTA

Step 3: Run Search Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run a basic search analysis (BLAST by default; use ``--search-tool mmseqs2`` for MMseqs2):

.. code-block:: bash

   eplace search query.fasta output_dir

This will:

1. Run BLAST search against the core_nt database
2. Filter results by identity (≥90%) and coverage (≥80%)
3. Extract representative sequences for each query at the genus level
4. Align sequences using MAFFT
5. Build phylogenetic trees using IQTree

View Results
~~~~~~~~~~~~

After the analysis completes, check the output directory:

.. code-block:: bash

   ls -R output_dir/

You'll find:

* ``blast_results.txt`` - Raw BLAST results
* ``blast_results_annotated.txt`` - BLAST results with taxonomic annotations
* One directory per query sequence containing:
  
  * Representative sequences (FASTA)
  * Multiple sequence alignment
  * Phylogenetic tree files

Common Workflows
----------------

Individual Analysis
~~~~~~~~~~~~~~~~~~~

Analyze each query sequence independently with custom parameters:

.. code-block:: bash

   eplace search query.fasta output_dir \
       --rank genus \
       --min-identity 95 \
       --min-coverage 85 \
       --num-threads 4

This creates one phylogenetic tree per query sequence.

Grouped Analysis
~~~~~~~~~~~~~~~~

Group queries by taxonomic classification:

.. code-block:: bash

   eplace grouped query.fasta output_dir \
       --rank genus \
       --group-rank family \
       --num-threads 4

This groups queries that match to the same family and creates one tree per group.

Relabel Existing Trees
~~~~~~~~~~~~~~~~~~~~~~~

Relabel an existing tree with taxonomic names at different ranks:

.. code-block:: bash

   # Relabel tree with genus names
   eplace relabel blast_results.txt input.treefile output_genus.treefile --rank genus
   
   # Relabel tree with species names (binomial nomenclature)
   eplace relabel blast_results.txt input.treefile output_species.treefile --rank species

This is useful when you want to create multiple versions of the same tree with different taxonomic labels without rebuilding the tree.

BLAST Only (No Alignment)
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you only want BLAST results without alignment and tree building:

.. code-block:: bash

   eplace search query.fasta output_dir --skip-alignment

High Stringency Search
~~~~~~~~~~~~~~~~~~~~~~

For more stringent matching:

.. code-block:: bash

   eplace search query.fasta output_dir \
       --min-identity 98 \
       --min-coverage 95

Custom Database Location
~~~~~~~~~~~~~~~~~~~~~~~~~

If your BLAST database is in a non-standard location:

.. code-block:: bash

   eplace search query.fasta output_dir \
       --blastdb-path /path/to/custom/blastdb

Using as a Python Library
--------------------------

You can also use ePLACE programmatically in your Python scripts.

Basic BLAST Workflow
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib import run_blast_search, process_blast_results_for_taxonomy

   # Run BLAST search with filtering
   success, filtered_hits = run_blast_search(
       query_fasta=Path("query.fasta"),
       output_file=Path("blast_results.txt"),
       min_identity=90.0,
       min_coverage=80.0
   )

   # Extract representative sequences by taxonomic rank
   results = process_blast_results_for_taxonomy(
       blast_hits=filtered_hits,
       output_dir=Path("output"),
       rank="genus"
   )

   # Print results
   for query_id, output_fasta in results.items():
       print(f"{query_id}: {output_fasta}")

Database Download
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from eplace_lib import setup_ncbi_database

   # Download the core_nt database
   success, message = setup_ncbi_database()
   print(f"Success: {success}, Message: {message}")

FASTA File Reading
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from eplace_lib.blast_analysis import FastaReader

   # Read sequences from FASTA file
   sequences = FastaReader.read_fasta(Path("input.fasta"))

   # Get sequence lengths
   lengths = FastaReader.get_sequence_lengths(Path("input.fasta"))

   for seq_id, length in lengths.items():
       print(f"{seq_id}: {length} bp")

Understanding Output Files
---------------------------

BLAST Results
~~~~~~~~~~~~~

* ``blast_results.txt`` - Tabular BLAST output with standard columns
* ``blast_results_annotated.txt`` - Same as above but with taxonomic annotations

Per-Query Directories
~~~~~~~~~~~~~~~~~~~~~

Each query sequence gets its own directory (``query_id/``) containing:

* ``query_id_representatives.fasta`` - Representative sequences selected by taxonomic rank
* ``query_id_with_query.fasta`` - Query sequence plus representatives
* ``query_id_trimmed.fasta`` - Sequences trimmed to aligned regions
* ``query_id_aligned.fasta`` - Multiple sequence alignment (MAFFT output)
* ``query_id_tree.treefile`` - Phylogenetic tree (Newick format)
* ``query_id_tree_labeled.treefile`` - Tree with taxonomic labels
* Additional IQTree output files

Grouped Analysis Output
~~~~~~~~~~~~~~~~~~~~~~~

For grouped analysis, you'll additionally see directories named by taxonomic group:

* ``Taxonomic_Group_Name/``
  
  * ``Taxonomic_Group_Name_combined.fasta`` - All queries and unique references
  * ``Taxonomic_Group_Name_trimmed.fasta`` - Trimmed sequences
  * ``Taxonomic_Group_Name_aligned.fasta`` - Multiple sequence alignment
  * ``Taxonomic_Group_Name_tree.treefile`` - Phylogenetic tree
  * Additional tree files

Choosing Between Workflows
---------------------------

Individual Workflow (``eplace search``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Use when:**

* You want to analyze each query in its own phylogenetic context
* Queries may be from diverse taxonomic groups
* You need separate trees for each sequence

**Advantages:**

* Independent analysis per query
* Clear interpretation per sequence
* No assumptions about relatedness

Grouped Workflow (``eplace grouped``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Use when:**

* You have multiple queries from related organisms
* You want to see queries together in phylogenetic context
* You want to reduce computational time for related sequences

**Advantages:**

* Combined phylogenetic analysis
* Fewer alignment/tree operations
* Better for comparative analysis

Command Reference
-----------------

For detailed command-line options, see:

* :doc:`cli` - Complete CLI reference
* ``eplace --help`` - General help
* ``eplace download --help`` - Download command help
* ``eplace search --help`` - Search command help
* ``eplace grouped --help`` - Grouped command help
* ``eplace relabel --help`` - Relabel command help

Next Steps
----------

* Read the detailed :doc:`workflows` documentation
* Learn about :doc:`blast_workflow` process
* Explore the :doc:`api` for programmatic access
* Check :doc:`ncbi_download` for database management

Troubleshooting
---------------

No BLAST hits found
~~~~~~~~~~~~~~~~~~~

If you get no BLAST hits:

1. Check that your sequences are in correct FASTA format
2. Try lowering ``--min-identity`` and ``--min-coverage`` thresholds
3. Verify your sequences are nucleotide sequences (not protein)
4. Ensure the BLAST database is properly installed

Command not found
~~~~~~~~~~~~~~~~~

If ``eplace`` command is not found:

1. Verify installation: ``pip show eplace``
2. Check that installation directory is in PATH
3. Try reinstalling: ``pip install --force-reinstall .``

Out of memory
~~~~~~~~~~~~~

If you run out of memory:

1. Process fewer sequences at a time
2. Use ``--skip-alignment`` to skip memory-intensive steps
3. Reduce ``--num-threads`` parameter
4. Consider using a machine with more RAM

Getting Help
------------

If you encounter issues:

1. Check this documentation
2. Review error messages carefully
3. Open an issue on GitHub: https://github.com/linsalrob/eplace/issues
4. Include error messages, command used, and system information
