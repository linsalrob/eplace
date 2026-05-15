# Installation

1. Create a conda environment for ePLACE:

```
mamba create -yn eplace 'python==3.12'
mamba activate eplace
```

2. Install the required dependencies

```
mamba install -y bioconda::blast bioconda::mmseqs2 bioconda::pytaxonkit bioconda::iqtree bioconda::mafft
```

> Note 1.
> At the time of writing there is an [issue](https://github.com/bioforensics/pytaxonkit/issues/50) with conda not installing the
> most current version of pytaxonkit if you are using python >=3.12. This code works with older versions of pytaxonkit.


> Note 2.
> You will need to download and set up the NCBI taxonomy databases for pytaxonkit; see the [taxonkit documentation](https://bioinf.shenwei.me/taxonkit/usage/#before-use) for detailed instructions of which NCBI taxonomy files to download.

3. Get and install eplace

```bash
pip install git+https://github.com/linsalrob/eplace.git
```

After installation, the `eplace` command will be available in your environment:

```bash
# Verify installation
eplace --help

# Show version
eplace --version
```

## Usage

Once installed, you can use the `eplace` command with four subcommands:

- `eplace download` - Download BLAST and/or MMseqs2 search databases
- `eplace search` - Run individual search workflow (one tree per query; BLAST or MMseqs2)
- `eplace grouped` - Run grouped search workflow (one tree per taxonomic group; BLAST or MMseqs2)
- `eplace relabel` - Relabel an existing tree with taxonomic names

For detailed help on each command:
```bash
eplace download --help
eplace search --help
eplace grouped --help
eplace relabel --help
```

Common download examples:
```bash
# BLAST core_nt only (default)
eplace download

# MMseqs2 NT only
eplace download --target mmseqs2

# MMseqs2 NT with taxonomy sidecars
eplace download --target mmseqs2 --add-taxonomy --ncbi-taxonomy /path/to/ncbi/taxonomy/current
```

See the [README.md](README.md) for complete documentation and examples.

