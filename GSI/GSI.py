#!/usr/bin/env python2

"""
This script calculates genealogical sorting index (GSI) for a given list of taxa. It iterates through all the trees in a file and outputs GSIs for each tree.

For more detail on GSI see: Cummings et al. 2008. A Genealogical Approach to Quantifying Lineage Divergence. Evolution, Vol. 62, No 9. pp.2411-2422

Files examples:

#tree.nwk:
((((a1,a2),(a3,a4)),((b1,b2),b3)),b4);
(((a1,a2),(a3,a4)),((b1,b2),(b3,b4)));
((((a1,a2),(b3,a4)),a3),((b1,b2),b4));
((((a1,a2),(a3,a4)),((b1,b2),b3)),b4);

#tree.names:
tree1
tree2
tree3
tree4

#GSI_tree.out:
TreeName    gr1 gr2
tree1   1.0 0.5625
tree2   1.0 1.0
tree3   0.5625  0.125
tree4   1.0 0.5625

#GSI_T_tree.out:
Group   GSI_T   p
gr1 0.890625    0.014000
gr2 0.562500    0.024000

#contact:
Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

# command:
python2 GSI.py -t tree.nwk -o tree.out -g "gr1[a1,a2,a3,a4];gr2[b1,b2,b3,b4]" -p 1000 -c tree.names

"""

############################# modules #############################

import argparse
import sys
import re
from ete2 import Tree
import matplotlib.pyplot as plt
import numpy as np
import random

############################# options #############################

class CommandLineParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = CommandLineParser()
parser.add_argument('-t', '--tree', help = 'file containing trees in newick format', type=str, required=True)
parser.add_argument('-c', '--coordinates', help = 'file indicating trees\' position on the genome', type=str, required=True)
parser.add_argument('-o', '--output', help = 'output name', type=str, required=True)
parser.add_argument('-g', '--groups', help = 'groups in the format "group1[taxon1,taxon2,taxon3];group2[taxon4,taxo5]"', type=str, required=True)
parser.add_argument('-p', '--permutation', help = 'number of permutations', type=int, required=False)
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
  GSI = (obsGS - minGS)/(maxGS-minGS) # formula 4 in Cummings et al. 2008
  '''
  maxGS = 1.0
  # count internal nodes
  dGr = count_subtree_internal_nodes(tree, group)
  n = len(group) - 1
  GS = float(n)/float(dGr)
  #print t.get_ascii(show_internal=True)  # for debugging
  #print n, dGr   # for debugging
  minGS = float(n)/float(count_all_internal_nodes(tree)) # modify this for large trees that include more leaves than specified in groups
  GSI = (GS - minGS) / (maxGS - minGS)
  return GSI

def tree_prop(tree, tree_file_name):
  ''' calculates the proportional representation of a tree in a nwk file.'''
  tsim = 0.0
  ttotal = 0.0
  treeFF = open(tree_file_name, 'r')
  for treeF in treeFF:
    tf = Tree(treeF)
    ttotal += 1.0
    # I use Robinson-Foulds metric to find the same trees.
    if Tree.compare(tree,tf)['norm_rf']== 0.0:
      tsim += 1.0
  treeFF.close()
  treeProb = (float(tsim/ttotal))
  return treeProb
  
############################# analysis #############################

# verify that tree and coordinates files are of the same length
if file_len(args.tree) != file_len(args.coordinates):
  raise IOError("%s and %s are of different length" % (args.tree, args.coordinates))

treeFile = open(args.tree, 'r')
coordinatesFile = open(args.coordinates, 'r')
outputFilet = open("GSI_"+args.output, 'w')
outputFileT = open("GSI_T_"+args.output, 'w')
groups = args.groups.split(';')

counter = 0

# create list of ID and names for each group
groupNames = []
indNames = []
GSIlists = []
gsiT = []
greter = []
treePlist = []
topologies = []

for i in range(len(groups)):
  groupsInd = re.split('\[|\]', groups[i])
  groupNames.append(groupsInd[0])
  indNames.append(groupsInd[1].split(","))
  GSIlists.append([]) # make lists to store GSI values for a histogram
  gsiT.append(0.0) # make lists to store GSI values for a gsiT
  greter.append(0) # make lists to store GSI values for the gsiT permutation list of lists

  
# make a header for the output
outputFilet.write("TreeName\t%s\n" % ('\t'.join(str(e) for e in groupNames)))

for tree, treeName in zip(treeFile,coordinatesFile):
  t = Tree(tree)
  #print t.get_ascii(show_internal=True)

  # verify that all names are present in a tree
  indNamesflat = [i for y in indNames for i in y]
  for indName in indNamesflat:
    if indName not in tree:
      raise ValueError("%s is not present in the tree %s" % (indName, tree))

  # resolve polytomies (probably unnecessary and can be disabled)
  t.resolve_polytomy(recursive=True)

  # add names to internal nodes for convenience of codding and debugging.
  assign_node_names(t)
  # print t.get_ascii(show_internal=True)  # for debugging

  # calculate GSI per tree
  GSIvalues = []
  for grn in range(len(groupNames)):
    GSIind = gsi(t, indNames[grn])
    GSIvalues.append(GSIind)
    GSIlists[grn].append(float(GSIind))

  # write output
  treeNameP = treeName.rstrip() # remove \n from a string
  GSIvaluesP = '\t'.join(str(e) for e in GSIvalues)
  outputFilet.write("%s\t%s\n" % (treeNameP, GSIvaluesP))

  # count topologies
  if not topologies:
    topologies.append(t)
    treeProb = tree_prop(t, args.tree)
    treePlist.append(treeProb)
  elif not any(Tree.compare(t,tt)['norm_rf'] == 0.0 for tt in topologies):
    topologies.append(t)
    treeProb = tree_prop(t, args.tree)
    treePlist.append(treeProb)

  # track the progress:
  counter += 1
  if counter % 100 == 0:
    print str(counter), "trees processed"

# calculate the GSI Total (formula 5 in Cummings et al. 2008)

# print topologies # for debugging
for top,p in zip(topologies, treePlist):
  for grnt in range(len(groupNames)):
    gsit = gsi(top, indNames[grnt])
    gsiT[grnt] += float(gsit)*p

############################# permutation of GSI Total #############################

outputFileT.write("Group\tGSI_T\tp\n") 
if args.permutation:
  pNperm = args.permutation
  for k in range(pNperm):
    gsiTp = []
    # permute labels
    for gr in groupNames:
      gsiTp.append(0.0) # make lists to store GSI values for the gsiT permutation

    indNamesFlat = [i for sl in indNames for i in sl]
    indNamesRand = random.sample(indNamesFlat, len(indNamesFlat))
    indNamesNew = []
    for l in range(len(indNames)):
      indNamesNew.append([])
      for ll in range(len(indNames[l])):
        indNamesNew[l].append(indNamesRand[ll])

    #calculate GSI with permuted labels
    for top,p in zip(topologies, treePlist):
      for grnt in range(len(groupNames)):
        gsitp = gsi(top, indNamesNew[grnt])
        gsiTp[grnt] += float(gsitp)*p

    # calculate p-values
    for i in range(len(gsiTp)):
      if gsiT[i] <= gsiTp[i]:
        greter[i] += 1

  # output GSI Total
  p = []
  for i in range(len(greter)):
    p = float(greter[i])/float(pNperm)
    outputFileT.write("%s\t%f\t%f\n" % (groupNames[i],gsiT[i], p))
else:
  for i in range(len(gsiT)):
    outputFileT.write("%s\t%f\tNA\n" % (groupNames[i], gsiT[i]))

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
outputFilet.close()
