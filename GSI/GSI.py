#!/usr/bin/env python2

"""
This script calculates genealogical sorting index (GSI) for a given list of taxa. It iterates through all the trees in a file and outputs GSIs for each tree.

For more detail on GSI see: Cummings et al. 2008. A Genealogical Approach to Quantifying Lineage Divergence. Evolution, Vol. 62, No 9. pp.2411-2422

Files examples:

#tree.nwk:
((((a1,a2),(a3,a4)),((b1,b2),b3)),b4);
(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));
((((a1,a2),(b3,a4)),a3),((b1,b2),b4));

#tree.names:
tree1
tree2
tree3

#tree.out:
TreeName	gr1	gr2
tree1	1.0	0.5625
tree2	1.0	1.0
tree3	0.5625	0.125


#contact:
Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

# command:
python2 GSI.py -t tree.nwk -o tree.out -g "gr1[a1,a2,a3,a4];gr2[b1,b2,b3,b4]" -c tree.names

"""

############################# modules #############################

import argparse
import sys
import re
from ete2 import Tree
import matplotlib.pyplot as plt
import numpy as np

############################# options #############################

class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-t', '--tree', help = 'file containing trees in newick format', type=str, required=True)
parser.add_argument('-c', '--coordinates', help = 'file indicating trees\' position on the genome', type=str, required=True)
parser.add_argument('-o', '--output', help = 'output name', type=str, required=True)
parser.add_argument('-g', '--groups', help = 'groups in the format "group1[taxon1,taxon2,taxon3];group2[taxon4,taxo5]"', type=str, required=True)
args = parser.parse_args()


############################# functions #############################

def file_len(filename):
  f = open(filename, 'r')
  for i, l in enumerate(f):
    pass
  f.close()
  return i + 1

def assign_node_names(tree):
  '''
  Assigns names to all internal nodes without names.
  The string "inter_" + number are used for names.
  '''
  myNodeName = 1
  for node in tree.traverse():
    if not node.name:
      node.add_features(name='inter_'+str(myNodeName))
      myNodeName += 1

def count_subtree_internal_nodes(tree, group):
  '''
  Counts internal nodes from leafs of group
  to the most recent common ancestor (MRCA) of the group
  '''
  ancestor = tree.get_common_ancestor(group).name # get the MRCA name
  nodes = [ancestor] # record the fist node 
  for group_name in group:
    leaf = tree.search_nodes(name=group_name)[0] # extract a leaf from a tree
    for a in leaf.iter_ancestors():
      if a.name != ancestor:
        nodes.append(a.name)
      else:
        break # break if the MRCA is reached
  #print set(nodes) # for debugging
  return int(len(set(nodes))) # count unique node names


def count_all_internal_nodes(tree):
  ''' Counts all internal nodes in a tree'''
  nodeL = 0
  for node in tree.traverse():
     if not node.is_leaf():
       nodeL +=1
  return nodeL

def gsi(tree, group):
  '''
  obsGS = n/sum(du-2)
  where d is the degree of node u of U total nodes uniting a group (estimated coalescent events) through the MRCA
  n is the minimum number of nodes (coalescent events) required to unite a group of size n + 1 through the MRCA
  maxGS is the maximum possible GS. The maxGS is reached when a group is monophyletic.
  minGS = n/sum(di-2)
  where i is one of I total nodes on the tree. Thus minGS would result if all nodes on a tree were required to unite a group.
  GSI = (obsGS - minGS)/(maxGS-minGS)
  '''
  maxGS = 1.0
  # count internal nodes
  dGr = count_subtree_internal_nodes(tree, group)
  n = len(group) - 1
  GS = float(n)/float(dGr)
  minGS = float(n)/float(count_all_internal_nodes(tree)) # modify this for large trees that include more leaves than specified in groups
  GSI = (GS - minGS) / (maxGS - minGS)
  return GSI

############################# analysis #############################

# verify that tree and coordinates files are of the same length
if file_len(args.tree) != file_len(args.coordinates):
  raise IOError("%s and %s are of different length" % (args.tree, args.coordinates))

treeFile = open(args.tree, 'r')
coordinatesFile = open(args.coordinates, 'r')
outputFile = open(args.output, 'w')
groups = args.groups.split(';')

counter = 0

# create list of ID and names for each group
groupNames = []
indNames = []
GSIlists = []
for i in range(len(groups)):
  groupsInd = re.split('\[|\]', groups[i])
  groupNames.append(groupsInd[0])
  indNames.append(groupsInd[1].split(","))
  GSIlists.append([]) # make lists to store GSI values for a histogram
#print groupNames # for debugging 
#print indNames # for debugging 

# make a header for the output
outputFile.write("TreeName\t%s\n" % ('\t'.join(str(e) for e in groupNames)))

for tree, treeName in zip(treeFile,coordinatesFile):
  t = Tree(tree)

  # verify that all names are present in a tree
  indNamesflat = [i for y in indNames for i in y]
  for indName in indNamesflat:
    if indName not in tree:
      raise ValueError("%s is not present in the tree %s" % (indName, tree))


  # add names to internal nodes for convenience of codding and debugging.
  assign_node_names(t)
  # print t.get_ascii(show_internal=True)  # for debugging

  # resolve polytomies
  t.resolve_polytomy(recursive=True)
  # print t.get_ascii(show_internal=True)  # for debugging
  # print '##############################\n' # for debugging

  # calculate GSI
  GSIvalues = []
  for i in range(len(groupNames)):
    GSIind = gsi(t, indNames[i])
    GSIvalues.append(GSIind)
    GSIlists[i].append(float(GSIind))

  # write output
  treeNameP = treeName.rstrip() # remove \n from a string
  GSIvaluesP = '\t'.join(str(e) for e in GSIvalues)
  outputFile.write("%s\t%s\n" % (treeNameP, GSIvaluesP))

  # track the progress:
  counter += 1
  if counter % 100 == 0:
    print str(counter), "trees processed"

############################# visualize #############################

#print GSI lists # for debugging

plt.hist(GSIlists, label=groupNames)
plt.xlim(0,1)
plt.xticks(np.arange(0,1.1,0.1))
plt.legend(loc = 0)
plt.title("GSI distribution", size = 20)
plt.xlabel("GSI")
plt.ylabel("Frequency")
plt.savefig(args.output +".png", dpi=90)

treeFile.close()
coordinatesFile.close()
outputFile.close()
