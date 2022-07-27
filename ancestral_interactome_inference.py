
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from ete3 import *


def check_output_dir(output_dir):
    output_dir = os.path.abspath(output_dir)
    # Create output directory if parent directory exists
    if not os.path.exists(output_dir):
        parent_dir = os.path.dirname(output_dir)
        if os.path.exists(parent_dir):
            os.mkdir(output_dir)
        else:
            sys.exit('EXIT: Unable to locate \'%s\'. Please provide a correct path to the output directory' % parent_dir)
        

# Fills:
# - an env. dictionary fams_to_edge that tells whether two input family names interact in at least one species and if such is case, returns the id of the edge 
# - an env. dictionary that associates to each edge_id, a list of mapped leaves (mapped taxa)
# - an env. dictionary that associates to each edge_id, a list of the two interacting families
# - an env. dictionary that tells whether an edge_id is present on a given leaf
# - a list of taxa/leaves considered in the analysis
# Returns:
# - an id_flag to start with when this function is called again for a novel PPI file
def fill_edges_data(filename, id_flag, leaves, fams_to_edge, edge_to_fams, edge_to_leaves, is_edge_on_leaf):
    fields = os.path.basename(filename).split('__')
    leaf = fields[0]
    leaves.append(leaf)
    with open(filename, mode = 'r') as f:
        for line in f:
            families = line.strip().split('\t')
            try:
                curr_id = fams_to_edge[families[0]][families[1]]
                # to deal with the duplicated edges
                if not leaf in edge_to_leaves[curr_id]:
                    edge_to_leaves[curr_id].append(leaf)
                    is_edge_on_leaf[curr_id][leaf] = True
            except:
                curr_id = id_flag
                if not families[0] in fams_to_edge: fams_to_edge[families[0]] = dict()
                if not families[1] in fams_to_edge: fams_to_edge[families[1]] = dict()
                if not families[0] in fams_to_edge[families[1]]: fams_to_edge[families[1]][families[0]] = dict()
                if not families[1] in fams_to_edge[families[0]]: fams_to_edge[families[0]][families[1]] = dict()
                fams_to_edge[families[0]][families[1]] = curr_id
                fams_to_edge[families[1]][families[0]] = curr_id
                edge_to_fams[curr_id] = families
                edge_to_leaves[curr_id] = list()
                is_edge_on_leaf[curr_id] = dict()
                id_flag += 1
                edge_to_leaves[curr_id].append(leaf)
                is_edge_on_leaf[curr_id][leaf] = True
    f.close()
    return id_flag, leaves


# This is just a wrapper of the above function
def load_edges(input_directory):
    id_flag = 0
    leaves = list()
    fams_to_edge = dict()
    edge_to_fams = dict()
    edge_to_leaves = dict()
    is_edge_on_leaf = dict()
    for filename in os.listdir(input_directory):
        id_flag, leaves = fill_edges_data(os.path.join(input_directory, filename), id_flag, leaves, 
                                          fams_to_edge, edge_to_fams, edge_to_leaves, is_edge_on_leaf)
    return leaves, edge_to_fams, edge_to_leaves, is_edge_on_leaf
        
    
# Name each internal node of a tree as following: <level_from_root>_<unique_identifier>
# Store additional node attributes such as nb of children and nb of attached leaves
def name_internal_nodes(node, node_id, nb_leaves, level_from_root):
    if node.is_leaf():
        nb_leaves = 1
        node.nb_children, node.nb_leaves, node.level = 0, 1, level_from_root
    else:
        node.name = str(level_from_root) + '_' + str(node_id)
        nb_leaves = 0
        node_id += 1
        children = node.get_children()
        for i in range(0, len(children)):
            flags = name_internal_nodes(children[i], node_id, nb_leaves, level_from_root + 1)
            node_id = flags[0]
            nb_leaves += flags[1]
        node.nb_children, node.nb_leaves, node.level = len(children), nb_leaves, level_from_root
    return node_id, nb_leaves
    

# Given a lineage and the edge_to_leaves mapping, 
# returns a dictionary that tells whether the edge is present or absent on any node of the lineage
def define_lost_nodes(node, fam, is_fam_on_leaf, losses, is_lost):
    if node.is_leaf():
        if node.name in is_fam_on_leaf[fam]:
            losses = 0
            is_lost[node.name] = False
        else:
            losses = 1
            is_lost[node.name] = True
    else:
        losses = 0
        children = node.get_children()
        for i in range(0, node.nb_children):
            losses += define_lost_nodes(children[i], fam, is_fam_on_leaf, losses, is_lost)
        if losses == node.nb_leaves:
            is_lost[node.name] = True
        else:
            is_lost[node.name] = False
    return losses


# Given a lineage and a dictionary that tells whether the edge is present or absent
# on any node of the lineage, fills:
# - a list of loss event(s) (node(s) from which the edge has been lost)
# - a list of all mapped nodes
def get_mapped_nodes(node, is_lost, loss_events, mapped_nodes):
    if is_lost[node.name]:
        loss_events.append(node.name)
    else:
        mapped_nodes.append(node.name)
        if not node.is_leaf():
            children = node.get_children()
            for i in range(0, node.nb_children):
                get_mapped_nodes(children[i], is_lost, loss_events, mapped_nodes)
                

# This function is just a wrapper of the two above functions
def map_edges_to_tree(output_file, edge_to_fams, edge_to_leaves, is_edge_on_leaf):
    edges = [*edge_to_leaves]
    mapping = dict()
    with open(output_file, mode = 'w+') as f:
        f.write("edge_id\tfamilyA\tfamilyB\tedge_starting_node\tedge_mapped_nodes\tedge_loss_events\n")
        for edge_id in edges:
            mapped_leaves = edge_to_leaves[edge_id]
            mapped_leaves.sort()
            key = ",".join(mapped_leaves)
            if key not in mapping:
                mapping[key] = dict()
                if len(mapped_leaves) == 1:
                    mapping[key]['starting_node_name'] = mapped_leaves[0]
                    mapping[key]['mapped_nodes'] = mapped_leaves[0]
                    mapping[key]['loss_events'] = ''
                else:
                    mapping[key]['starting_node'] = tree.get_common_ancestor(mapped_leaves)
                    mapping[key]['starting_node_name']  = mapping[key]['starting_node'].name
                    is_lost = dict()
                    define_lost_nodes(mapping[key]['starting_node'], edge_id, is_edge_on_leaf, 0, is_lost)
                    loss_events = list()
                    mapped_nodes = list()
                    get_mapped_nodes(mapping[key]['starting_node'], is_lost, loss_events, mapped_nodes)
                    mapping[key]['mapped_nodes'] = ','.join(mapped_nodes)
                    mapping[key]['loss_events'] = ','.join(loss_events)
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % 
                    (edge_id,
                     edge_to_fams[edge_id][0],
                     edge_to_fams[edge_id][1],
                     mapping[key]['starting_node_name'],
                     mapping[key]['mapped_nodes'],
                     mapping[key]['loss_events']))
    f.close()
            
        
        
# Read arguments
pkg_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
parser = argparse.ArgumentParser(description='Map edges of PPIs to nodes of a reference phylogenetic tree')
parser.add_argument('-t', dest='tree', type=str, help='path to the reference phylogenetic tree in Newick')
parser.add_argument('-e', dest='edges_dir', type=str, help='path to directory where networks/edges are stored, each under the name <taxon_name_on_the_tree>__<whatever_suffix>')
parser.add_argument('-o', dest='outdir', type = str, help='path to the output dir')
args = parser.parse_args()
    
# LOAD INPUT DATA
print('loading input data ...')
check_output_dir(args.outdir)
leaves, edge_to_fams, edge_to_leaves, is_edge_on_leaf = load_edges(args.edges_dir)
tree = Tree(args.tree, format = 1)

# NAME AND ANNOTATE INTERNAL NODES
print('naming internal nodes ...')
name_internal_nodes(tree, 0, 0, 0)
tree.write(features=['name'], format=1, outfile=os.path.join(args.outdir, "tree_with_internal_nodes.txt"))

# MAP EACH EDGE ON THE TREE
print('definining all mapped nodes + all loss event(s) of each edge ...')
map_edges_to_tree(os.path.join(args.outdir, 'datation.tsv'), edge_to_fams, edge_to_leaves, is_edge_on_leaf)
