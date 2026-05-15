Command-Line Interface
======================

The ``eplace`` command provides a unified interface to all ePLACE functionality through three subcommands.

Installation Verification
--------------------------

After installing the package with ``pip install .`` or ``pip install -e .``, the ``eplace`` command will be available:

.. code-block:: bash

   # Verify installation
   eplace --help

   # Check version
   eplace --version

Commands Overview
-----------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Command
     - Description
   * - ``eplace download``
     - Download BLAST and/or MMseqs2 search databases
   * - ``eplace search``
     - Run individual search workflow (one tree per query; BLAST by default, MMseqs2 via ``--search-tool``)
   * - ``eplace grouped``
     - Run grouped BLAST workflow (one tree per taxonomic group)
   * - ``eplace relabel``
     - Relabel phylogenetic tree with taxonomic names

eplace download
---------------

Download and setup the NCBI core_nt BLAST database and/or MMseqs2 NT database.

Usage
~~~~~

.. code-block:: bash

   eplace download [--target {blast,mmseqs2,both}] [--force] [MMSEQS_OPTIONS]

Options
~~~~~~~

.. option:: --target {blast,mmseqs2,both}

   Database backend(s) to download

   Default: ``blast``

.. option:: --force

   Force redownload even if database exists

.. option:: --mmseqs-db-dir PATH

   Path to MMseqs2 database root directory. Defaults to ``$MMSEQS_DB_DIR``,
   then ``$MMSEQS2DB``, then ``~/mmseqs2db``.

.. option:: --mmseqs-threads INT

   Number of threads for MMseqs2 download/taxonomy commands

.. option:: --add-taxonomy

   Add taxonomy sidecar files to MMseqs2 NT database after download.

.. option:: --ncbi-taxonomy PATH

   Path to NCBI taxonomy dump directory containing ``nodes.dmp``,
   ``names.dmp``, and ``merged.dmp``. Required with ``--add-taxonomy``.

.. option:: --acc2taxid-dir PATH

   Path to accession2taxid files. Defaults to ``$ACC2TAXID_DIR`` or
   ``<ncbi-taxonomy>/accession2taxid``.

.. option:: --taxonomy-workdir PATH

   Working directory for MMseqs taxonomy mapping files.

.. option:: --skip-memory-check

   Skip RAM preflight checks for MMseqs2 download/taxonomy.

Examples
~~~~~~~~

.. code-block:: bash

   # Download BLAST database to default location ($BLASTDB or ~/blastdb)
   eplace download

   # Download MMseqs2 NT database
   eplace download --target mmseqs2 --mmseqs-db-dir /path/to/mmseqs_db

   # Download MMseqs2 NT database and add taxonomy sidecar files
   eplace download --target mmseqs2 --add-taxonomy --ncbi-taxonomy /path/to/ncbi/taxonomy/current

Notes
~~~~~

* BLAST DB location: ``$BLASTDB`` or ``~/blastdb``
* MMseqs2 DB location: ``$MMSEQS_DB_DIR``, then ``$MMSEQS2DB``, or ``~/mmseqs2db``
* MMseqs2 NT download typically requires at least 64 GiB RAM
* MMseqs2 taxonomy integration typically requires at least 128 GiB RAM

eplace search
------------

Run sequence search with individual taxonomy analysis. Creates one phylogenetic tree per query sequence.

Usage
~~~~~

.. code-block:: bash

   eplace search QUERY_FASTA OUTPUT_DIR [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

.. option:: QUERY_FASTA

   Path to query FASTA file containing sequences to search

.. option:: OUTPUT_DIR

   Output directory for results (will be created if it doesn't exist)

Optional Arguments
~~~~~~~~~~~~~~~~~~

Taxonomy Options
^^^^^^^^^^^^^^^^

.. option:: --rank {phylum,class,order,family,genus,species}

   Taxonomic rank for representative selection
   
   Default: ``genus``

.. option:: --tree-label-rank {phylum,class,order,family,genus,species}

   Taxonomic rank for tree labeling
   
   Default: ``genus``

Filtering Options
^^^^^^^^^^^^^^^^^

.. option:: --min-identity FLOAT

   Minimum percent identity for BLAST hits
   
   Default: ``90.0``

.. option:: --min-coverage FLOAT

   Minimum query coverage percentage
   
   Default: ``80.0``

Database Options
^^^^^^^^^^^^^^^^

.. option:: --database NAME

   BLAST database name
   
   Default: ``core_nt``

.. option:: --blastdb-path PATH

   Path to BLAST database directory

Search Backend Options
^^^^^^^^^^^^^^^^^^^^^^

.. option:: --search-tool {blast,mmseqs2}

   Sequence search tool to use

   Default: ``blast``

.. option:: --blast-db-source LABEL

   Provenance label written to ``search_metadata.json`` for reproducibility.
   If not provided, this defaults to ``--database``.

.. option:: --mmseqs-database NAME

   MMseqs2 database name. Defaults to ``--database`` when omitted.
   Only used with ``--search-tool mmseqs2``.

.. option:: --mmseqs-db-path PATH

   Path to MMseqs2 database directory. Only used with
   ``--search-tool mmseqs2``. Defaults to ``$MMSEQS_DB_DIR``, then
   ``$MMSEQS2DB``, then ``~/mmseqs2db``.

.. option:: --mmseqs-db-source LABEL

   Provenance label for MMseqs2 database origin, written to
   ``search_metadata.json``. Only used with ``--search-tool mmseqs2``.

.. option:: --mmseqs-sensitivity FLOAT

   MMseqs2 sensitivity (range 1 to 7.5).
   Lower values are faster but less sensitive; higher values are slower but
   more thorough.
   Typical tuning: use around ``3-4`` for faster exploratory runs, and
   ``6-7`` when you need more comprehensive hit recovery.
   Only used with ``--search-tool mmseqs2``.

   Default: ``5.7``

.. option:: --mmseqs-search-type INT

   MMseqs2 ``--search-type`` value passed to ``mmseqs easy-search``.
   Common values include ``2`` (translated), ``3`` (nucleotide), and
   ``4`` (translated nucleotide backtrace).
   Only used with ``--search-tool mmseqs2``.

   Default: ``3``

Performance Options
^^^^^^^^^^^^^^^^^^^

.. option:: --num-threads INT

   Number of threads for BLAST and alignment
   
   Default: ``1``

Workflow Options
^^^^^^^^^^^^^^^^

.. option:: --overwrite-existing-blast

   Overwrite existing BLAST results

.. option:: --skip-alignment

   Skip alignment and tree building steps

.. option:: --output-classification PATH

   Path to output classification TSV file

Examples
~~~~~~~~

.. code-block:: bash

   # Basic usage with default parameters
   eplace search query.fasta output_dir

   # With custom parameters
   eplace search query.fasta output_dir \
       --rank genus \
       --min-identity 95 \
       --min-coverage 85 \
       --num-threads 4

   # Skip alignment and tree building (BLAST only)
   eplace search query.fasta output_dir --skip-alignment

   # Use custom BLAST database location
   eplace search query.fasta output_dir --blastdb-path /path/to/blastdb

   # Use MMseqs2 with explicit backend configuration
   eplace search query.fasta output_dir \
       --search-tool mmseqs2 \
       --mmseqs-db-path /path/to/mmseqsdb \
       --mmseqs-database NT

Output Structure
~~~~~~~~~~~~~~~~

.. code-block:: text

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
   └── ...

eplace grouped
--------------

Run sequence search with grouped taxonomy analysis. Groups queries by taxonomic rank and creates one phylogenetic tree per group.

Usage
~~~~~

.. code-block:: bash

   eplace grouped QUERY_FASTA OUTPUT_DIR [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

.. option:: QUERY_FASTA

   Path to query FASTA file containing sequences to search

.. option:: OUTPUT_DIR

   Output directory for results (will be created if it doesn't exist)

Optional Arguments
~~~~~~~~~~~~~~~~~~

Taxonomy Options
^^^^^^^^^^^^^^^^

.. option:: --rank {phylum,class,order,family,genus,species}

   Taxonomic rank for representative selection
   
   Default: ``genus``

.. option:: --group-rank {phylum,class,order,family,genus,species}

   Taxonomic rank for grouping sequences
   
   Default: ``class``

.. option:: --tree-label-rank {phylum,class,order,family,genus,species}

   Taxonomic rank for tree labeling
   
   Default: ``genus``

.. option:: --combined-tree-label-rank {phylum,class,order,family,genus,species}

   Taxonomic rank for labeling the combined tree (optional)
   
   Default: Not set (combined tree will not be built)
   
   The grouped workflow can create a combined tree from all groups, but this is optional because it can be very time-consuming with large datasets. If you want to build the combined tree, specify this parameter with the desired taxonomic rank for labeling. If not provided, only individual group trees will be built.

Filtering Options
^^^^^^^^^^^^^^^^^

.. option:: --min-identity FLOAT

   Minimum percent identity for BLAST hits
   
   Default: ``90.0``

.. option:: --min-coverage FLOAT

   Minimum query coverage percentage
   
   Default: ``80.0``

Database Options
^^^^^^^^^^^^^^^^

.. option:: --database NAME

   BLAST database name
   
   Default: ``core_nt``

.. option:: --blastdb-path PATH

   Path to BLAST database directory

Search Backend Options
^^^^^^^^^^^^^^^^^^^^^^

.. option:: --search-tool {blast,mmseqs2}

   Sequence search tool to use

   Default: ``blast``

.. option:: --blast-db-source LABEL

   Provenance label written to ``search_metadata.json`` for reproducibility.
   If not provided, this defaults to ``--database``.

.. option:: --mmseqs-database NAME

   MMseqs2 database name. Defaults to ``--database`` when omitted.
   Only used with ``--search-tool mmseqs2``.

.. option:: --mmseqs-db-path PATH

   Path to MMseqs2 database directory. Only used with
   ``--search-tool mmseqs2``. Defaults to ``$MMSEQS_DB_DIR``, then
   ``$MMSEQS2DB``, then ``~/mmseqs2db``.

.. option:: --mmseqs-db-source LABEL

   Provenance label for MMseqs2 database origin, written to
   ``search_metadata.json``. Only used with ``--search-tool mmseqs2``.

.. option:: --mmseqs-sensitivity FLOAT

   MMseqs2 sensitivity (range 1 to 7.5).
   Lower values are faster but less sensitive; higher values are slower but
   more thorough.
   Typical tuning: use around ``3-4`` for faster exploratory runs, and
   ``6-7`` when you need more comprehensive hit recovery.
   Only used with ``--search-tool mmseqs2``.

   Default: ``5.7``

.. option:: --mmseqs-search-type INT

   MMseqs2 ``--search-type`` value passed to ``mmseqs easy-search``.
   Common values include ``2`` (translated), ``3`` (nucleotide), and
   ``4`` (translated nucleotide backtrace).
   Only used with ``--search-tool mmseqs2``.

   Default: ``3``

Performance Options
^^^^^^^^^^^^^^^^^^^

.. option:: --num-threads INT

   Number of threads for BLAST and alignment
   
   Default: ``1``

Workflow Options
^^^^^^^^^^^^^^^^

.. option:: --overwrite-existing-blast

   Overwrite existing BLAST results

.. option:: --skip-alignment

   Skip alignment and tree building steps

.. option:: --alignment-tolerance INT

   Maximum coordinate difference for alignment consistency
   
   Default: ``50``

.. option:: --output-classification PATH

   Path to output classification TSV file

Examples
~~~~~~~~

.. code-block:: bash

   # Basic usage (group by class, default)
   eplace grouped query.fasta output_dir

   # Group by different taxonomic rank
   eplace grouped query.fasta output_dir --group-rank order

   # Specify both representative and grouping ranks
   eplace grouped query.fasta output_dir --rank genus --group-rank family

   # Skip alignment and tree building
   eplace grouped query.fasta output_dir --skip-alignment

   # Use MMseqs2 backend for grouped workflow
   eplace grouped query.fasta output_dir \
       --search-tool mmseqs2 \
       --mmseqs-db-path /path/to/mmseqsdb \
       --mmseqs-database NT

Output Structure
~~~~~~~~~~~~~~~~

.. code-block:: text

   output_dir/
   ├── blast_results.txt              # Raw BLAST results
   ├── blast_results_annotated.txt    # BLAST results with taxonomic annotations
   ├── query1_id/                     # Per-query representative sequences
   │   └── query1_id_representatives.fasta
   ├── Taxonomic_Group_1/             # One directory per taxonomic group
   │   ├── Taxonomic_Group_1_combined.fasta        # All queries + unique references
   │   ├── Taxonomic_Group_1_trimmed.fasta         # Trimmed to aligned regions
   │   ├── Taxonomic_Group_1_aligned.fasta         # Multiple sequence alignment
   │   ├── Taxonomic_Group_1_tree.treefile         # Phylogenetic tree
   │   ├── Taxonomic_Group_1_tree_labeled.treefile # Tree with taxonomic labels
   │   └── Taxonomic_Group_1_tree.* (other IQTree files)
   ├── combined_all_groups_trimmed.fasta           # Combined alignment of all groups
   ├── combined_all_groups_aligned.fasta           # Multiple sequence alignment
   ├── combined_all_groups_tree.treefile           # Combined phylogenetic tree
   ├── combined_all_groups_tree_labeled.treefile   # Combined tree with taxonomic labels
   └── ...

eplace relabel
--------------

Relabel a phylogenetic tree with taxonomic names from BLAST results. This command allows you to replace sequence IDs in an existing tree with taxonomic names at any specified rank.

Usage
~~~~~

.. code-block:: bash

   eplace relabel BLAST_OUTPUT TREE_FILE OUTPUT_TREE [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

.. option:: BLAST_OUTPUT

   Path to BLAST output file (tabular format with taxonomy)
   
   The BLAST results file should contain taxonomic information for the sequences in the tree.

.. option:: TREE_FILE

   Path to input tree file (Newick format)
   
   The phylogenetic tree to be relabeled with taxonomic names.

.. option:: OUTPUT_TREE

   Path to output relabeled tree file
   
   The new tree file with taxonomic labels will be written to this path.

Optional Arguments
~~~~~~~~~~~~~~~~~~

Taxonomy Options
^^^^^^^^^^^^^^^^

.. option:: --rank {phylum,class,order,family,genus,species}

   Taxonomic rank for tree labeling
   
   Default: ``genus``
   
   When using ``species``, the output will use binomial nomenclature (genus + species).

Database Options
^^^^^^^^^^^^^^^^

.. option:: --blastdb-path PATH

   Path to BLAST database directory
   
   Optional parameter, not required for relabeling operation.

Examples
~~~~~~~~

.. code-block:: bash

   # Relabel tree with genus names (default)
   eplace relabel blast_results.txt input.treefile output_labeled.treefile
   
   # Relabel tree with species names (genus + species binomial)
   eplace relabel blast_results.txt input.treefile output_species.treefile --rank species
   
   # Relabel tree with family names
   eplace relabel blast_results.txt input.treefile output_family.treefile --rank family
   
   # Relabel tree with order names
   eplace relabel blast_results.txt input.treefile output_order.treefile --rank order

Key Features
~~~~~~~~~~~~

* **Flexible Taxonomic Ranks**: Supports all standard taxonomic ranks from phylum to species
* **Binomial Nomenclature**: Species rank uses "genus species" format for proper scientific names
* **Topology Preservation**: Maintains the original tree structure while updating labels
* **Format Compatibility**: Works with standard Newick format trees
* **Reversed Sequences**: Handles sequences with _R_ prefix (from MAFFT orientation)
* **Label Cleaning**: Automatically cleans labels for Newick format compatibility

Use Cases
~~~~~~~~~

The ``eplace relabel`` command is useful in several scenarios:

1. **Re-labeling at Different Ranks**
   
   Generate multiple versions of the same tree with different taxonomic granularity without rebuilding the tree:
   
   .. code-block:: bash
   
      # Create genus-level tree
      eplace relabel blast_results.txt tree.treefile tree_genus.treefile --rank genus
      
      # Create family-level tree from same data
      eplace relabel blast_results.txt tree.treefile tree_family.treefile --rank family

2. **External Tree Tools**
   
   Add taxonomic labels to trees generated by external phylogenetic tools:
   
   .. code-block:: bash
   
      # After building tree with RAxML, IQTree, or FastTree
      eplace relabel blast_results.txt external_tree.nwk labeled_tree.nwk --rank genus

3. **Visualization Preparation**
   
   Create publication-ready trees with appropriate taxonomic labels:
   
   .. code-block:: bash
   
      # For species-level visualization
      eplace relabel blast_results.txt tree.treefile publication_tree.treefile --rank species

4. **Updating Taxonomy**
   
   Update tree labels when taxonomy databases are updated or corrected:
   
   .. code-block:: bash
   
      # Re-run BLAST with updated database, then relabel
      eplace relabel new_blast_results.txt old_tree.treefile updated_tree.treefile

Notes
~~~~~

* The BLAST results file must contain taxonomic information for the sequences
* Tree file must be in Newick format (standard phylogenetic tree format)
* Sequence IDs in the tree must match accession numbers in the BLAST results
* Labels are automatically cleaned to comply with Newick format requirements (spaces, colons, parentheses are replaced)
* If taxonomy information is missing for a sequence, it will be skipped with a warning

Output
~~~~~~

The output is a Newick format tree file with sequence IDs replaced by taxonomic names at the specified rank.

.. code-block:: text

   Input tree:
   ((NZ_CP012345:0.01,NZ_CP067890:0.02):0.03,(NC_012345:0.01,NC_067890:0.02):0.03);
   
   Output tree (--rank genus):
   ((Escherichia:0.01,Salmonella:0.02):0.03,(Bacillus:0.01,Staphylococcus:0.02):0.03);
   
   Output tree (--rank species):
   ((Escherichia_coli:0.01,Salmonella_enterica:0.02):0.03,(Bacillus_subtilis:0.01,Staphylococcus_aureus:0.02):0.03);

Workflow Comparison
-------------------

Individual Workflow (eplace search)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Best for:** Analyzing each query sequence in its own phylogenetic context

**Process:**

1. Run BLAST search for all queries
2. Extract representative sequences for each query at specified rank
3. Create one directory per query
4. Build one alignment and tree per query

**Output:** Separate phylogenetic trees for each query sequence

Grouped Workflow (eplace grouped)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Best for:** Analyzing multiple related queries together in a single phylogenetic context

**Process:**

1. Run BLAST search for all queries
2. Extract representative sequences for each query
3. Group queries by specified taxonomic rank (e.g., class, order)
4. Combine all queries in a group with unique reference sequences
5. Build one alignment and tree per group

**Output:** Phylogenetic trees with multiple queries grouped by taxonomy

Common Use Cases
----------------

Quick BLAST search without trees
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   eplace search query.fasta results --skip-alignment

Relabel existing trees
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Relabel tree with different taxonomic ranks
   eplace relabel blast_results.txt input.treefile genus_tree.treefile --rank genus
   eplace relabel blast_results.txt input.treefile family_tree.treefile --rank family
   
   # Use with trees from external tools
   eplace relabel blast_results.txt raxml_tree.nwk labeled_tree.nwk --rank species

High stringency search
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   eplace search query.fasta results \
       --min-identity 95 \
       --min-coverage 90

Multi-threaded analysis
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   eplace search query.fasta results --num-threads 8

Group related sequences
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   eplace grouped query.fasta results \
       --rank species \
       --group-rank genus

Troubleshooting
---------------

Command not found
~~~~~~~~~~~~~~~~~

If ``eplace`` command is not found after installation:

.. code-block:: bash

   # Check if it's installed
   pip show eplace

   # Reinstall
   pip install --force-reinstall .

   # Or add to PATH
   export PATH="$HOME/.local/bin:$PATH"

Dependencies missing
~~~~~~~~~~~~~~~~~~~~

Some features require external tools:

**BLAST+**: Required for all workflows

.. code-block:: bash

   sudo apt-get install ncbi-blast+  # Ubuntu/Debian
   brew install blast                 # macOS

**TaxonKit**: Required for taxonomy lookups

.. code-block:: bash

   conda install -c bioconda taxonkit

**MAFFT**: Required for alignment (unless --skip-alignment)

.. code-block:: bash

   sudo apt-get install mafft  # Ubuntu/Debian
   brew install mafft          # macOS

**IQTree**: Required for tree building (unless --skip-alignment)

.. code-block:: bash

   sudo apt-get install iqtree  # Ubuntu/Debian
   brew install iqtree          # macOS

See Also
--------

* :doc:`installation` - Installation instructions
* :doc:`workflows` - Detailed workflow documentation
* :doc:`blast_workflow` - Complete workflow guide
* :doc:`ncbi_download` - Database management details
