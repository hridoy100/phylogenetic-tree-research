# import dendropy
from io import StringIO
import dendropy
from dendropy import TaxonSet, Tree, TreeList, Node
from dendropy.simulate import treesim
import pickle
from collections import OrderedDict

def getAllTrees(filename):
    # taxa = dendropy.TaxonNamespace(["a", "b", "c", "d"])
    trees = TreeList.get_from_path(filename, schema="newick")
    return trees

def getSingleTree(filename):
    # taxa = dendropy.TaxonNamespace(["a", "b", "c", "d"])
    tree = Tree.get_from_path(filename, schema="newick")
    return tree


# trees = getAllTrees("sample.input")
trees = getAllTrees("output_paup.greedy.tree")

clustersWithCount = dict()
f = open("clusters_paup.txt", "w")
# f = open("clusters.txt", "w")
f.write(str(len(trees))+'\n')
allIntNodes = []
for tree in trees:
    intNodes = tree.internal_nodes()
    f.write(str(len(intNodes)-1) + '\n')

    # print(tree.leaf_nodes())
    allTaxon = ""

    for idx, node in enumerate(tree.leaf_nodes()):
        splitStrings = str(node.taxon).split("'")
        actualTaxonLabel = splitStrings[1]
        if idx>0:
            allTaxon = allTaxon+" "+actualTaxonLabel
        else:
            allTaxon = allTaxon + actualTaxonLabel

    # print("all Leaves ", allTaxon)
    # print("all Internal Nodes ", intNodes)
    internalNodesExcludingRoot = []
    for index, node in enumerate(intNodes):
        # print("node: ", node)
        leaves = node.leaf_nodes()
        # print("leaves ", leaves)
        label = ""
        for idx, leaf in enumerate(leaves):
            splitStrings = str(leaf.taxon).split("'")
            actualTaxonLabel = splitStrings[1]  #instead of storing 'a' we extract only a

            if idx>0:
                label = label+" "+actualTaxonLabel
            else:
                label = label+actualTaxonLabel
            # print(label)
            # print(leaf.taxon)
        # print("label: ", label)
        node.label = label
        # print("node label: ", node.description())
        if node.label != allTaxon:
            f.write(str(node.label)+'\n')
            node.label = node.label.replace(" ", "")
            node.label = "".join(sorted(node.label))
            internalNodesExcludingRoot.append(node)

    allIntNodes.append(internalNodesExcludingRoot)
# f.write(allIntNodes)

f.close()

