# Ancestral interactome inference

## What it does

For each edge present in at least two extant networks (on the basis of (ortho|homo)logy)

```
1: present / 0: absent

          -- 1 X--O (edge present in the network of this species)         
       --|
      |   -- 1 X--O (edge present in the network of this species)
    --|
   |  |   -- 0
 --|   --|
   |      -- 0 
   |
    -------- 1 X--O (edge present in the network of this species)
```

returns the inferred evolutionary history of the edge along the tree (using parsimony)
This allows ancestral networks inferences. 

```
          -- 1 X--O
      1--|
      |   -- 1 X--O
   1--|
   |  |   -- 0
1--|  0--|
   |      -- 0
   |
    -------- 1 X--O
```

## Installation

```bash
pip install ete3
cd <your_install_directory>
clone https://github.com/TeamAIRE/ancestral_interactome_inference.git
```

## Usage

Run on linux, with python3

```
usage: ancestral_interactome_inference.py [-h] [-t TREE] [-e EDGES_DIR]
                                          [-o OUTDIR]

Map edges of PPIs to nodes of a reference phylogenetic tree

optional arguments:
  -h, --help    show this help message and exit
  -t TREE       path to the reference phylogenetic tree in Newick
  -e EDGES_DIR  path to directory where networks/edges are stored, each under
                the name <taxon_name_on_the_tree>__<whatever_suffix>
  -o OUTDIR     path to the output dir
```

## Format of input data

see ```ancestral_interactome_inference/first_lines_of_input_data_files_to_visualize_format``` for template

* TREE

A phylogenetic tree of the taxa for which interactome data is available (in newick)

* EDGES_DIR

A directory where the interactomes are stored.

Each interactome file must be named according to the corresponding taxon name (as referenced on the aformentionned tree) followed by "__" and whatever suffix.
An interactome file is a tabular file (tab-separated tsv file) with no header and two columns (nodeA, nodeB). Each line corresponds to an edge. Each node (protein) must be named according to the orthologous group id in order to allow comparisons of interactomes between extant taxa. 