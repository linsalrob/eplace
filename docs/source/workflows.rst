Workflows
=========

ePLACE provides two main workflow types for analyzing environmental DNA sequences: **Individual** and **Grouped** workflows.

Workflow Overview
-----------------

ePLACE provides two main workflow types for analyzing environmental DNA sequences: **Individual** and **Grouped** workflows. Additionally, the **Relabel** command allows you to post-process phylogenetic trees with taxonomic labels.

Common Pipeline Steps
~~~~~~~~~~~~~~~~~~~~~

1. **Sequence Search** - Search query sequences with BLAST or MMseqs2
2. **Filter Results** - Apply identity and coverage thresholds
3. **Taxonomic Classification** - Classify hits using TaxonKit
4. **Representative Selection** - Select representative sequences per taxonomic rank
5. **Sequence Extraction** - Extract sequences from BLAST database
6. **Sequence Trimming** - Trim to aligned regions based on BLAST coordinates
7. **Multiple Sequence Alignment** - Align using MAFFT (optional)
8. **Phylogenetic Tree Building** - Build trees using IQTree (optional)
9. **Tree Relabeling** - Relabel trees with taxonomic names (optional or standalone)

Individual Workflow
-------------------

The individual workflow (``eplace search``) processes each query sequence independently.

When to Use
~~~~~~~~~~~

* Analyzing diverse sequences that may not be related
* Need independent phylogenetic context for each sequence
* Want clear, per-sequence results
* Queries span multiple taxonomic groups

Process
~~~~~~~

.. code-block:: text

   For each query sequence:
   1. Sequence search (BLAST/MMseqs2) → filter hits
   2. Classify hits taxonomically
   3. Select representatives at specified rank
   4. Create query-specific directory
   5. Extract and trim sequences
   6. Align sequences (query + representatives)
   7. Build phylogenetic tree

Output Structure
~~~~~~~~~~~~~~~~

.. code-block:: text

   output_dir/
   ├── blast_results.txt
   ├── blast_results_annotated.txt
   ├── query1/
   │   ├── query1_representatives.fasta
   │   ├── query1_with_query.fasta
   │   ├── query1_trimmed.fasta
   │   ├── query1_aligned.fasta
   │   └── query1_tree.treefile
   ├── query2/
   │   └── ...
   └── ...

Usage Example
~~~~~~~~~~~~~

.. code-block:: bash

   # Basic individual workflow
   eplace search queries.fasta output_dir

   # With custom parameters
   eplace search queries.fasta output_dir \
       --rank genus \
       --min-identity 95 \
       --min-coverage 85 \
       --num-threads 4

   # Search only, no alignment/tree building
   eplace search queries.fasta output_dir --skip-alignment

Grouped Workflow
----------------

The grouped workflow (``eplace grouped``) combines queries by taxonomic classification for joint analysis.

When to Use
~~~~~~~~~~~

* Analyzing related sequences from similar taxonomic groups
* Want to see multiple queries in single phylogenetic context
* Need comparative analysis of related sequences
* Want to reduce computational overhead

Process
~~~~~~~

.. code-block:: text

   1. Sequence search all queries (BLAST/MMseqs2) → filter hits
   2. Classify hits taxonomically
   3. Select representatives for each query at specified rank
   4. Group queries by taxonomic rank (e.g., class, order, family)
   5. For each group:
      a. Combine all queries in group
      b. Collect unique reference sequences
      c. Remove redundant references
      d. Verify alignment consistency
      e. Extract and trim sequences
      f. Align all sequences together
      g. Build single phylogenetic tree
   6. Build combined tree from all groups (optional):
      a. Combine representatives from all groups
      b. Align combined sequences
      c. Build phylogenetic tree showing relationships across all groups

Output Structure
~~~~~~~~~~~~~~~~

.. code-block:: text

   output_dir/
   ├── blast_results.txt
   ├── blast_results_annotated.txt
   ├── query1/
   │   └── query1_representatives.fasta
   ├── query2/
   │   └── query2_representatives.fasta
   ├── Bacteria_Proteobacteria/
   │   ├── Bacteria_Proteobacteria_combined.fasta
   │   ├── Bacteria_Proteobacteria_trimmed.fasta
   │   ├── Bacteria_Proteobacteria_aligned.fasta
   │   └── Bacteria_Proteobacteria_tree.treefile
   ├── Bacteria_Firmicutes/
   │   └── ...
   ├── combined_all_groups_trimmed.fasta           # Combined from all groups
   ├── combined_all_groups_aligned.fasta           # Alignment of combined sequences
   ├── combined_all_groups_tree.treefile           # Combined phylogenetic tree
   ├── combined_all_groups_tree_labeled.treefile   # Combined tree with labels
   └── ...

Usage Example
~~~~~~~~~~~~~

.. code-block:: bash

   # Basic grouped workflow (groups by class)
   eplace grouped queries.fasta output_dir

   # Group by different rank
   eplace grouped queries.fasta output_dir --group-rank order

   # Specify both representative and grouping ranks
   eplace grouped queries.fasta output_dir \
       --rank genus \
       --group-rank family

   # With custom parameters
   eplace grouped queries.fasta output_dir \
       --rank genus \
       --group-rank family \
       --min-identity 95 \
       --min-coverage 85 \
       --num-threads 4

Relabel Workflow
----------------

The relabel workflow (``eplace relabel``) allows you to post-process existing phylogenetic trees by replacing sequence IDs with taxonomic names.

When to Use
~~~~~~~~~~~

* You have an existing phylogenetic tree that needs taxonomic labels
* You want to create multiple versions of a tree with different taxonomic ranks
* You built a tree with external tools (RAxML, FastTree, etc.) and want to add taxonomic labels
* You need to update tree labels after taxonomy database updates
* You want to avoid rebuilding trees just to change label granularity

Process
~~~~~~~

.. code-block:: text

   1. Parse BLAST results to get taxonomy information
   2. Read existing phylogenetic tree (Newick format)
   3. Create mapping of sequence IDs to taxonomic names
   4. Replace sequence IDs with taxonomic labels at specified rank
   5. Clean labels for Newick format compatibility
   6. Write relabeled tree to output file

Required Inputs
~~~~~~~~~~~~~~~

* **BLAST Results**: File containing BLAST hits with taxonomy information
* **Tree File**: Existing phylogenetic tree in Newick format
* **Output Path**: Where to save the relabeled tree

Output
~~~~~~

A new Newick format tree file with taxonomic labels replacing sequence accession numbers.

Usage Examples
~~~~~~~~~~~~~~

.. code-block:: bash

   # Basic relabeling with genus names (default)
   eplace relabel blast_results.txt input.treefile output.treefile
   
   # Relabel with species names (binomial nomenclature)
   eplace relabel blast_results.txt input.treefile species_tree.treefile --rank species
   
   # Create multiple versions at different ranks
   eplace relabel blast_results.txt tree.treefile genus_tree.treefile --rank genus
   eplace relabel blast_results.txt tree.treefile family_tree.treefile --rank family
   eplace relabel blast_results.txt tree.treefile order_tree.treefile --rank order
   
   # Use with trees from external tools
   eplace relabel blast_results.txt raxml_tree.nwk labeled_tree.nwk --rank genus

Advantages
~~~~~~~~~~

* **Fast**: No tree rebuilding required
* **Flexible**: Easily create trees with different taxonomic granularity
* **Compatible**: Works with trees from any phylogenetic tool
* **Preserves Topology**: Maintains original tree structure
* **Format Support**: Handles standard Newick format and MAFFT-oriented sequences

Comparison
----------

+-----------------------+-------------------------+-------------------------+-------------------------+
| Feature               | Individual              | Grouped                 | Relabel                 |
+=======================+=========================+=========================+=========================+
| Analysis unit         | Per query               | Per taxonomic group     | Existing tree           |
+-----------------------+-------------------------+-------------------------+-------------------------+
| Trees generated       | One per query           | One per group           | N/A (modifies existing) |
+-----------------------+-------------------------+-------------------------+-------------------------+
| Alignment scope       | Query + its refs        | All queries + unique    | N/A                     |
|                       |                         | refs in group           |                         |
+-----------------------+-------------------------+-------------------------+-------------------------+
| Best for              | Diverse sequences       | Related sequences       | Post-processing trees   |
+-----------------------+-------------------------+-------------------------+-------------------------+
| Computational cost    | Higher (more trees)     | Lower (fewer trees)     | Minimal (no rebuild)    |
+-----------------------+-------------------------+-------------------------+-------------------------+
| Interpretation        | Independent per query   | Comparative across      | Visual/taxonomic labels |
|                       |                         | queries                 |                         |
+-----------------------+-------------------------+-------------------------+-------------------------+
| BLAST required        | Yes                     | Yes                     | Yes (for taxonomy)      |
+-----------------------+-------------------------+-------------------------+-------------------------+
| Tree building         | Yes                     | Yes                     | No (uses existing)      |
+-----------------------+-------------------------+-------------------------+-------------------------+

Workflow Parameters
-------------------

Taxonomic Rank Selection
~~~~~~~~~~~~~~~~~~~~~~~~~

Both workflows use taxonomic ranks to organize results:

* ``--rank``: Level at which to select representative sequences

  * Options: ``phylum``, ``class``, ``order``, ``family``, ``genus``, ``species``
  * Default: ``genus``
  * Higher ranks = fewer, more diverse representatives
  * Lower ranks = more, more specific representatives

* ``--group-rank`` (grouped only): Level at which to group queries

  * Options: ``phylum``, ``class``, ``order``, ``family``, ``genus``, ``species``
  * Default: ``class``
  * Determines how queries are combined
  * Should typically be equal to or higher than ``--rank``

* ``--tree-label-rank``: Level at which to label tree tips

  * Options: ``phylum``, ``class``, ``order``, ``family``, ``genus``, ``species``
  * Default: ``genus``
  * Determines taxonomic labels on phylogenetic tree

* ``--combined-tree-label-rank`` (grouped only): Level at which to label the combined tree (optional)

  * Options: ``phylum``, ``class``, ``order``, ``family``, ``genus``, ``species``
  * Default: Not set (combined tree will not be built)
  * Controls labels on the combined tree built from all groups
  * Can be different from ``--tree-label-rank`` for more flexibility
  * **Note:** Building the combined tree can be very time-consuming with large datasets, so it is only built when explicitly requested

Filtering Parameters
~~~~~~~~~~~~~~~~~~~~

Control which BLAST hits are included:

* ``--min-identity``: Minimum percent identity (default: 90.0)

  * Range: 0-100
  * Higher = more stringent matching
  * Typical values: 90-99

* ``--min-coverage``: Minimum query coverage (default: 80.0)

  * Range: 0-100
  * Percentage of query sequence covered by alignment
  * Higher = require longer alignments

Database Parameters
~~~~~~~~~~~~~~~~~~~

Configure search database usage:

* ``--database``: BLAST database name (default: core_nt)

  * Standard databases: nt, core_nt, refseq_rna, etc.
  * Must be installed in BLASTDB path

* ``--blastdb-path``: Custom path to BLAST databases

  * Overrides $BLASTDB environment variable
  * Use for non-standard locations

* ``--search-tool``: Search backend to use (default: ``blast``)

  * Options: ``blast``, ``mmseqs2``
  * Applies to individual and grouped workflows
  * Choose MMseqs2 for alternate search backend behavior

* ``--mmseqs-database``: MMseqs2 database name

  * Defaults to ``--database`` when omitted
  * Used only when ``--search-tool mmseqs2`` is selected

* ``--mmseqs-db-path``: Custom path to MMseqs2 database directory

  * Defaults to ``$MMSEQS_DB_DIR``, then ``$MMSEQS2DB``, then ``~/mmseqs2db``
  * Used only when ``--search-tool mmseqs2`` is selected

* ``--blast-db-source`` / ``--mmseqs-db-source``: Provenance labels for search metadata

  * Recorded in ``search_metadata.json`` for reproducibility
  * Use these labels to document custom or non-standard databases

Performance Parameters
~~~~~~~~~~~~~~~~~~~~~~

Optimize analysis speed:

* ``--num-threads``: Number of CPU threads to use (default: 1)

  * Used by BLAST, MAFFT, and IQTree
  * Set to number of available cores for speed
  * Higher values = faster but more resource intensive

Workflow Control
~~~~~~~~~~~~~~~~

Fine-tune workflow behavior:

* ``--skip-alignment``: Skip MAFFT alignment and IQTree tree building

  * Useful for BLAST-only analysis
  * Saves time and computational resources
  * Still performs extraction and trimming

* ``--overwrite-existing-blast``: Re-run BLAST even if results exist

  * By default, existing BLAST results are reused
  * Use when BLAST parameters change

* ``--alignment-tolerance`` (grouped only): Maximum coordinate difference for alignment consistency (default: 50)

  * Checks if BLAST hits align to similar regions
  * Larger values = more permissive grouping
  * Smaller values = stricter alignment requirements

* ``--output-classification``: Path to output classification TSV file

  * Generates summary table of taxonomic classifications
  * Useful for downstream analysis

Advanced Usage
--------------

Combining Workflows with Relabel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use relabel to create multiple tree versions from workflow outputs:

.. code-block:: bash

   # Run workflow once
   eplace search queries.fasta output --rank genus
   
   # Create trees with different label ranks
   cd output/query1
   eplace relabel ../../output/blast_results_annotated.txt query1_tree.treefile query1_genus.treefile --rank genus
   eplace relabel ../../output/blast_results_annotated.txt query1_tree.treefile query1_family.treefile --rank family
   eplace relabel ../../output/blast_results_annotated.txt query1_tree.treefile query1_species.treefile --rank species

This approach is efficient because you only build the tree once, then relabel it multiple times.

Chaining Analyses
~~~~~~~~~~~~~~~~~

You can run different analyses on the same BLAST results:

.. code-block:: bash

   # First run: search only (no alignment/tree building)
   eplace search queries.fasta output1 --skip-alignment

   # Analyze at different ranks using same BLAST results
   eplace search queries.fasta output2 --rank genus
   eplace search queries.fasta output3 --rank species

Batch Processing
~~~~~~~~~~~~~~~~

Process multiple query files:

.. code-block:: bash

   for file in queries/*.fasta; do
       base=$(basename "$file" .fasta)
       eplace search "$file" "output_$base" --num-threads 8
   done

Custom Filtering Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~

Apply different stringency levels:

.. code-block:: bash

   # High stringency
   eplace search queries.fasta output_strict \
       --min-identity 98 \
       --min-coverage 95

   # Medium stringency
   eplace search queries.fasta output_medium \
       --min-identity 95 \
       --min-coverage 85

   # Low stringency
   eplace search queries.fasta output_relaxed \
       --min-identity 90 \
       --min-coverage 75

Hierarchical Grouping
~~~~~~~~~~~~~~~~~~~~~~

Analyze at multiple taxonomic levels:

.. code-block:: bash

   # Group by class
   eplace grouped queries.fasta output_class --group-rank class

   # Group by order (more specific)
   eplace grouped queries.fasta output_order --group-rank order

   # Group by family (even more specific)
   eplace grouped queries.fasta output_family --group-rank family

Best Practices
--------------

Choosing Parameters
~~~~~~~~~~~~~~~~~~~

1. **Start with defaults** to get a baseline
2. **Adjust identity/coverage** based on sequence quality
3. **Choose rank based on analysis goals**:

   * Phylum/Class: Broad overview
   * Order/Family: Balanced detail
   * Genus/Species: Specific identification

4. **Use grouped workflow** when queries are known to be related
5. **Use ``--num-threads``** to speed up analysis

Quality Control
~~~~~~~~~~~~~~~

1. Check BLAST results for reasonable number of hits
2. Verify taxonomic classifications make sense
3. Inspect alignments for quality
4. Review phylogenetic trees for expected topology

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

1. Use ``--skip-alignment`` for exploratory analysis
2. Increase ``--num-threads`` on multi-core systems
3. Use grouped workflow to reduce number of trees
4. Filter query sequences before analysis
5. Consider splitting very large query files

Troubleshooting
---------------

No sequences in group
~~~~~~~~~~~~~~~~~~~~~

If grouped workflow creates empty groups:

* Queries may be too diverse
* Try lower ``--group-rank`` (e.g., phylum instead of class)
* Check taxonomic classifications of queries

Alignment consistency errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If grouped workflow reports alignment inconsistencies:

* BLAST hits align to different regions of references
* Increase ``--alignment-tolerance``
* Or use individual workflow instead

Tree building fails
~~~~~~~~~~~~~~~~~~~

If IQTree fails to build tree:

* Check alignment quality
* Ensure sufficient sequences (≥3)
* Verify sequences have overlapping regions
* Try ``--skip-alignment`` to see raw sequences

See Also
--------

* :doc:`cli` - Complete command-line reference
* :doc:`blast_workflow` - Detailed BLAST workflow documentation
* :doc:`api` - Python API for custom workflows
