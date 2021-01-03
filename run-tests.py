import time

import numpy as np
from dendropy import TaxonSet, Tree, TreeList, Node
import random
from datetime import datetime

DEBUG = False
PRINT_ESSENTIALS = False


def get_all_trees(filename):
    tree_list = TreeList.get_from_path(filename, schema="newick")
    return tree_list


def get_single_tree(filename):
    single_tree = Tree.get_from_path(filename, schema="newick")
    return single_tree


def is_compatible(c1, c2):
    if DEBUG:
        print('intersection: ', c1, ', ', c2, ' => \n',
              c1.intersection(c2), 'length: ', len(c1.intersection(c2)))
    # check pairwise compatibility
    if not (c1.issubset(c2)) and not (c1.issuperset(c2)) and len(
            c1.intersection(c2)) > 0:
        return False
    return True


def get_pairwise_table(cluster_list, cluster_count):
    num_clusters = len(cluster_list)
    table = np.zeros((num_clusters, num_clusters))
    for i in range(0, num_clusters - 1):
        for j in range(i + 1, num_clusters):
            c1 = cluster_list[i]
            c2 = cluster_list[j]
            if is_compatible(c1, c2):
                table[i][j] += 1 * cluster_count[c1] * cluster_count[c2]  # or other values
            else:
                table[i][j] += -1  # or other values

    return np.sum(table[table >= 0])


def get_subsets(frozen_set, cluster_count):
    subsets_of_a_set = []
    for sets in cluster_count:
        if frozen_set.issuperset(sets):
            subsets_of_a_set.append(sets)
    return subsets_of_a_set


def intelligent_tree_reconstruct(cluster_file):
    start_time = time.time()
    f = open(cluster_file, "r")
    n_trees = int(f.readline())
    cluster_count = dict()
    # Counting clusters--------------------------------------------------------
    for _ in range(n_trees):
        # for each tree input number of cluster
        n_cluster = int(f.readline())
        # takes the number of clusters
        for _ in range(n_cluster):
            cluster_list = f.readline().split()
            cluster = frozenset(cluster_list)
            cluster_count[cluster] = cluster_count.get(cluster, 0) + 1
    unsorted_cluster_count = cluster_count
    cluster_count = dict(sorted(cluster_count.items(), key=lambda kv: kv[1], reverse=True))
    duplicates = {}
    for k, v in enumerate(cluster_count):
        duplicates.setdefault(cluster_count[v], []).append(v)

    if DEBUG:
        print("duplicate clusters according to counts: ", duplicates)
        print("cluster_count", cluster_count)

    if PRINT_ESSENTIALS:
        print("duplicate clusters according to counts: ", duplicates)

    # y is our desired set of clusters ------------------------------------------------------------------------------
    y = dict()

    # For each duplicate clusters, make table of compatibility and store it....
    tables = dict()
    for cluster in cluster_count:
        subsets = get_subsets(cluster, cluster_count)
        tables[cluster] = get_pairwise_table(subsets, cluster_count)

    cluster_count_all = cluster_count.copy()

    for iter, c in enumerate(cluster_count.keys()):
        if DEBUG:
            print('iteration---> ', iter, ' duplicate clusters for count: (', cluster_count[c], ")  => \n",
                  duplicates.get(cluster_count[c]))

        if PRINT_ESSENTIALS:
            print('iteration (', iter, ') ----------------------------------------------------------------')

        # here, subsetCount is a dictionary which counts the number of subsets a single set has in the cluster_count
        # dictionary without choosing a set randomly, we count the subsets under each duplicate count sets.
        # then we choose the set with maximum count
        subset_count = dict()
        # for duplicate in duplicates.get(c[1]):
        for duplicate in duplicates.get(cluster_count[c]):
            subset_count[duplicate] = tables[duplicate]

        if DEBUG:
            print("subset count: ", subset_count)
        subset_count = sorted(subset_count.items(), key=lambda kv: kv[1], reverse=False)
        if DEBUG:
            print("count of subsets: ", subset_count)
        # The last item has the maximum subset count. so, that item is popped.
        # again, if the subset counts are found equal, we have no way but to choose randomly........
        max_subset = subset_count.pop()
        if DEBUG:
            print("chosen subset frozenSet with max count: ", max_subset[0])

        if PRINT_ESSENTIALS:
            print("chosen subset frozenSet with max count: ", max_subset)

        # we need to remove the item 'maxsubsetcount' from duplicates dictionary-------------------------------
        duplicates[cluster_count[c]] = []
        for subsets in subset_count:
            duplicates[cluster_count[c]].append(subsets[0])

        if DEBUG:
            print("new duplicates for count: (", cluster_count[c], ") => ", duplicates)

        if PRINT_ESSENTIALS:
            print("new duplicates for count: (", cluster_count[c], ") => ", duplicates)

        if len(y) == 0:
            y[max_subset[0]] = cluster_count[c]
            if PRINT_ESSENTIALS:
                print("selected subset-------------- ", max_subset)

            del cluster_count_all[max_subset[0]]
        else:
            # For all clusters in y check if c is pairwise compatible with all of them
            flag = 0
            for cPrime in y:
                if not is_compatible(max_subset[0], cPrime):
                    flag = 1
                    break
            # if it is compatible then add into the final cluster set--------------------------------------------
            if flag == 0:
                y[max_subset[0]] = cluster_count[c]
                if PRINT_ESSENTIALS:
                    print("selected subset-------------- ", max_subset)
                del cluster_count_all[max_subset[0]]

        if PRINT_ESSENTIALS:
            print('----------------------------------------------------------------\n')

    return y, start_time


def random_selection(cluster_file):
    start_time = time.time()
    f = open(cluster_file, "r")
    n_trees = int(f.readline())
    cluster_count = dict()
    # Counting clusters--------------------------------------------------------
    for _ in range(n_trees):
        # for each tree input number of cluster
        n_cluster = int(f.readline())
        # takes the number of clusters
        for _ in range(n_cluster):
            cluster_list = f.readline().split()
            cluster = frozenset(cluster_list)
            cluster_count[cluster] = cluster_count.get(cluster, 0) + 1

    cluster_count = sorted(cluster_count.items(), key=lambda kv: kv[1], reverse=True)
    # duplicates is a dictionary which has counts as keys and clusters as values--------------------------------------
    duplicates = {}
    for k, v in enumerate(cluster_count):
        duplicates.setdefault(v[1], []).append(v[0])

    if DEBUG:
        print("duplicate clusters according to counts: ", duplicates)
        print("cluster_count", cluster_count)

    # y is our desired set of clusters -----------------------------------------------------
    y = dict()
    for c in cluster_count:
        if DEBUG:
            print('duplicate clusters for count: (', c[1], ")  => \n", duplicates.get(c[1]))

        duplicates_for_this_count = duplicates.get(c[1])
        # random.seed(0)
        rand_index = random.randint(0, len(duplicates_for_this_count) - 1)

        if DEBUG:
            print("random index: ", rand_index)

        selected_random_duplicate_for_this_count = duplicates_for_this_count[rand_index]

        del duplicates.get(c[1])[rand_index]
        if DEBUG:
            print("selected Duplicate: ", selected_random_duplicate_for_this_count)

        if len(y) == 0:
            y[selected_random_duplicate_for_this_count] = c[1]
        else:
            # For all clusters in y check if c is pairwise compatible with all of them
            flag = 0
            for cPrime in y:
                # print('isSubSet: ', c[0], ', ', cPrime, c[0].issubset(cPrime))
                if DEBUG:
                    print('intersection: ', selected_random_duplicate_for_this_count, ', ', cPrime, ' => \n',
                          selected_random_duplicate_for_this_count.intersection(cPrime), 'length: ',
                          len(selected_random_duplicate_for_this_count.intersection(cPrime)))
                # check pairwise compatibility
                if not (selected_random_duplicate_for_this_count.issubset(cPrime)) and not (
                        selected_random_duplicate_for_this_count.issuperset(cPrime)) and len(
                    selected_random_duplicate_for_this_count.intersection(cPrime)) > 0:
                    flag = 1
                    break
            # if it is compatible then add into the final cluster set------------------------------------
            if flag == 0:
                y[selected_random_duplicate_for_this_count] = c[1]

    return y, start_time


def normal_tree_reconstruct(cluster_file):
    start_time = time.time()
    f = open(cluster_file, "r")
    n_trees = int(f.readline())
    cluster_count = dict()
    # Counting clusters--------------------------------------------------------
    for _ in range(n_trees):
        # for each tree input number of cluster
        n_cluster = int(f.readline())

        # takes the number of clusters
        for _ in range(n_cluster):
            cluster_list = f.readline().split()
            cluster = frozenset(cluster_list)
            cluster_count[cluster] = cluster_count.get(cluster, 0) + 1
            # the get function return 0 if cluster does not exists already

    cluster_count = dict(sorted(cluster_count.items(), key=lambda kv: kv[1], reverse=True))
    # duplicates is a dictionary which has counts as keys and clusters as values------------------------
    duplicates = {}
    for k, v in enumerate(cluster_count):
        duplicates.setdefault(cluster_count[v], []).append(v)

    if DEBUG:
        print("duplicate clusters according to counts: ", duplicates)
        print("cluster_count", cluster_count)

    # y is our desired set of clusters ----------------------------------------------------------
    y = dict()
    for c in cluster_count:
        if DEBUG:
            print('duplicate clusters for count: (', cluster_count[c], ")  => \n", duplicates.get(cluster_count[c]))

        # here, subsetCount is a dictionary which counts the number of subsets a single set has in the cluster_count
        # dictionary without choosing a set randomly, we count the subsets under each duplicate count sets.
        # then we choose the set with maximum count
        subset_count = dict()
        for duplicate in duplicates.get((cluster_count[c])):
            num_of_subsets = 0
            frozen_set = duplicate
            for sets in cluster_count:
                # if frozenSet.issuperset(sets[0]):
                if frozen_set.issuperset(sets):
                    num_of_subsets += 1
            subset_count[duplicate] = num_of_subsets

        # subsetCount is sorted in increasing order------------------------------
        subset_count = sorted(subset_count.items(), key=lambda kv: kv[1], reverse=False)
        if DEBUG:
            print("countOfSubsets: ", dict(subset_count), type(subset_count))
        # The last item has the maximum subset count. so, that item is popped.
        # again, if the subset counts are found equal, we have no way but to choose randomly........
        # max_subset = subset_count.pop()
        # max_counted_set = subset_count.pop()
        max_subset = subset_count[-1]
        same_count_keys = [k for k, v in dict(subset_count).items() if v == max_subset[1]]
        if len(same_count_keys) > 1:
            rand_index = random.randint(0, len(same_count_keys) - 1)
            max_subset = subset_count[rand_index]
            subset_count.pop(rand_index)
            if DEBUG:
                print("chosen random index: ", rand_index)
        else:
            subset_count.pop()
        if DEBUG:
            print("Keys with same count: ", same_count_keys)
            print("chosen subset frozenSet with max count: ", max_subset[0])

        # we need to remove the item 'max_subset_count' from duplicates dictionary-------------------------------
        duplicates[cluster_count[c]] = []
        for subsets in subset_count:
            duplicates[cluster_count[c]].append(subsets[0])

        if DEBUG:
            print("new duplicates for count: (", cluster_count[c], ") => ", duplicates)

        if len(y) == 0:
            y[max_subset[0]] = cluster_count[c]
        else:
            # For all clusters in y check if c is pairwise compatible with all of them
            flag = 0
            for cPrime in y:
                if not is_compatible(max_subset[0], cPrime):
                    flag = 1
                    break
            # if it is compatible then add into the final cluster set-------------------------
            if flag == 0:
                y[max_subset[0]] = cluster_count[c]

    return y, start_time

"""
# input_file = open("drive/My Drive/Thesis/tree-input1.input", "r")
input_file = open("tree-input1.input", "r")
# input_file = open("sample50.input", "r")
input_rows = {}
for index, line in enumerate(input_file):
    input_rows[index] = line  # strip removes \n
if DEBUG:
    print(input_rows)
input_file.close()
"""
# input_file_name = "test-inputs.input"

# input_file_name = "DS.input"
cluster_file_name = "clusters.txt"

n_random_selection = 100
n_tree_reconstruct = 1
# output_result_file = "output50trees.txt"
now = datetime.now()
current_time = now.strftime("%H-%M-%S")

# for cnt in ['50', '100']:
for cnt in ['100']:
    output_result_file = "Output/output " + cnt + " trees " + current_time + ".txt"
    output_file = open(output_result_file, "w")
    print("Running for ", cnt, " input trees")
    output_file.write("Running for " + cnt + " input trees\n")
    # for index in range(0, len(input_rows)):
    # for dataset in range(1, 6):
    for dataset in range(4, 5):
        # for index in range(0, 5):
        # print("removed input tree # ", index+1)
        # output_file.write("removed input tree # " + str(index+1) + "\n")
        # file = open(input_file_name, "w")
        # file = open("drive/My Drive/Thesis/Datasets/"+cnt+"_input/DS"+cnt+"_"+dataset, "w")
        # for key, value in input_rows.items():
        # #     if key != index:
        #       file.write(value)
        # file.close()
        trees = get_all_trees("Datasets/" + cnt + "_input/DS" + cnt + "_" + str(dataset) + ".txt")
        clustersWithCount = dict()
        clusterFile = open(cluster_file_name, "w")
        clusterFile.write(str(len(trees)) + '\n')
        allIntNodes = []
        for tree in trees:
            intNodes = tree.internal_nodes()
            clusterFile.write(str(len(intNodes) - 1) + '\n')
            allTaxon = ""
            for idx, node in enumerate(tree.leaf_nodes()):
                splitStrings = str(node.taxon).split("'")
                actualTaxonLabel = splitStrings[1]
                if idx > 0:
                    allTaxon = allTaxon + " " + actualTaxonLabel
                else:
                    allTaxon = allTaxon + actualTaxonLabel
            internalNodesExcludingRoot = []
            for node in intNodes:
                leaves = node.leaf_nodes()
                label = ""
                for leaf_idx, leaf in enumerate(leaves):
                    splitStrings = str(leaf.taxon).split("'")
                    actualTaxonLabel = splitStrings[1]  # instead of storing 'a' we extract only a
                    if leaf_idx > 0:
                        label = label + " " + actualTaxonLabel
                    else:
                        label = label + actualTaxonLabel
                node.label = label
                if node.label != allTaxon:
                    clusterFile.write(str(node.label) + '\n')
                    node.label = node.label.replace(" ", "")
                    node.label = "".join(sorted(node.label))
                    internalNodesExcludingRoot.append(node)

            allIntNodes.append(internalNodesExcludingRoot)

        clusterFile.close()

        print("--------------------------------------------------------")
        output_file.write("--------------------------------------------------------"+"\n")
        print("Iteration: " + str(dataset))
        output_file.write("Iteration: " + str(dataset) + "\n")
        print("Intelligent Tree Reconstruct =====>")
        output_file.write("Intelligent Tree Reconstruct =====>" + "\n")

        max_cluster_length = 0
        max_cluster = {}
        cluster_length_probability = {}
        start_time = time.time()
        for i in range(n_tree_reconstruct):
            # y, _ = normal_tree_reconstruct(cluster_file_name)
            y, _ = intelligent_tree_reconstruct(cluster_file_name)
            if len(y) > max_cluster_length:
                max_cluster = y
                max_cluster_length = len(y)
            cluster_length_probability[len(y)] = cluster_length_probability.get(len(y), 0) + 1
        for key, value in cluster_length_probability.items():
            cluster_length_probability[key] = (value / n_tree_reconstruct) * 100

        cluster_length_probability = dict(sorted(cluster_length_probability.items(), reverse=True))

        print("final cluster: ", max_cluster)
        print("final cluster length: ", max_cluster_length)
        print("cluster probabilities: ", cluster_length_probability)
        print("--- %s seconds ---" % (time.time() - start_time))
        # print("--------------------------------------------------------")

        output_file.write("final cluster: " + str(max_cluster) + "\n")
        output_file.write("final cluster length: " + str(max_cluster_length) + "\n")
        output_file.write("cluster probabilities: " + str(cluster_length_probability) + "\n")
        output_file.write("--- %s seconds ---" % (time.time() - start_time) + "\n")

        # output_file.write("--------------------------------------------------------")

        # y, start_time = intelligent_tree_reconstruct(cluster_file_name)
        # y, start_time = normal_tree_reconstruct(cluster_file_name)
        # print("final cluster: ", y)
        # print("final cluster length: ", len(y))
        # print("--- %s seconds ---" % (time.time() - start_time))
        # print("Random Selection:")
        #
        # output_file.write("final cluster: " + str(y) + "\n")
        # output_file.write("final cluster length: " + str(len(y)) + "\n")
        # output_file.write("--- %s seconds ---" % (time.time() - start_time) + "\n")
        print("\n" + "Random Selection =====>")
        output_file.write("\n" + "Random Selection =====>" + "\n")

        max_cluster_length = 0
        expectation = 0
        max_cluster = {}
        cluster_length_probability = {}
        start_time = time.time()
        for i in range(n_random_selection):
            y, _ = random_selection(cluster_file_name)

            """
            if len(y) > max_cluster_length:
                max_cluster = y
                max_cluster_length = len(y)
            cluster_length_probability[len(y)] = cluster_length_probability.get(len(y), 0) + 1
            """
            cluster_length_probability[len(y)] = cluster_length_probability.get(len(y), 0) + 1
            if len(y) * cluster_length_probability[len(y)] > expectation:
                max_cluster = y
                expectation = len(y) * cluster_length_probability[len(y)]
                max_cluster_length = len(y)

        for key, value in cluster_length_probability.items():
            cluster_length_probability[key] = (value / n_random_selection) * 100

        cluster_length_probability = dict(sorted(cluster_length_probability.items(), reverse=True))

        print("final cluster: ", max_cluster)
        print("final cluster length: ", max_cluster_length)
        print("cluster probabilities: ", cluster_length_probability)
        print("--- %s seconds ---" % (time.time() - start_time))
        print("--------------------------------------------------------\n")

        output_file.write("final cluster: " + str(max_cluster) + "\n")
        output_file.write("final cluster length: " + str(max_cluster_length) + "\n")
        output_file.write("cluster probabilities: " + str(cluster_length_probability) + "\n")
        output_file.write("--- %s seconds ---" % (time.time() - start_time) + "\n")
        output_file.write("--------------------------------------------------------\n")
        # time.sleep(50)
output_file.close()
