#!/usr/bin/env python
import copy

class Lineage():
    def __init__(self, time, left_child, right_child):
        # time since the start of the simulation when this lineage coalesced
        self.time = time
        self.left = left_child
        self.right = right_child
        # sparse array storing mutation events that occured to this lineage
        self.mutations = {}

    def __str__(self):
        return "Lineage coalescing at " + str(self.time)

    def mutate(self, n_mutations):
        self.mutations[n_mutations] = 1

def save(tree, filehandle):
    pickle.dump(tree, filehandle)

def load(filehandle):
    return pickle.load(filehandle)

def branch_length(lineage):
    """Calculate the total branch length of the tree rooted at lineage."""
    if lineage.left is None or lineage.right is None:
        # this is an initial lineage
        return 0
    else:
        left = lineage.left
        right = lineage.right
        t_left = (lineage.time - left.time) + branch_length(left)
        t_right = (lineage.time - right.time) + branch_length(right)
        return t_left + t_right

def get_leaves(root):
    """Return all the leaves under a root"""
    if root.left == None and root.right == None:
        return [root]
    elif root.left == None:
        return get_leaves(root.right)
    elif root.right == None:
        return get_leaves(root.left)
    else:
        leaves = get_leaves(root.left)
        leaves.extend(get_leaves(root.right))
        return leaves

def pairwise_times(root):
    """Calculate T_ij for tree rooted at root for i < j"""
    if root.left == None and root.right == None:
        return []
    elif root.left == None:
        return pairwise_times(root.right)
    elif root.right == None:
        return pairwise_times(root.left)
    else:
        lefts = len(get_leaves(root.left))
        rights = len(get_leaves(root.right))
        # For each left-right pair, add this coalescence time.
        # with n left leaves and m right leaves, there are nm pairs that go
        # through this node.
        times = [root.time for i in range(lefts * rights)]
        times.extend(pairwise_times(root.left))
        times.extend(pairwise_times(root.right))
        return times

def build_sequences_rec(lineage, mutations):
    m = copy.copy(mutations)
    m.update(lineage.mutations)
    if lineage.left is not None and lineage.right is not None:
        lefts = build_sequences_rec(lineage.left, m)
        rights = build_sequences_rec(lineage.right, m)
        lefts.extend(rights)
        return lefts
    else:
        # we're a leaf node (or the tree is broken)
        return [m]

def build_sequences(lineage):
    return build_sequences_rec(lineage, {})

def convert_sequence(mutations, n_mutations):
    sequence = []
    for i in range(0, n_mutations):
        if i in mutations:
            sequence.append(1)
        else:
            sequence.append(0)
    return sequence

def build_site_frequency(sequences):
    n = len(sequences)
    n_mutations = len(sequences[0])
    frequencies = {}
    def add(count):
        if count in frequencies:
            frequencies[count] += 1
        else:
            frequencies[count] = 1
    for i in range(0, n_mutations):
        count = 0
        for sequence in sequences:
            count += sequence[i]
            # adding the minimum gives us the folded frequency spectrum
            add(min(count, n-count))
    return frequencies
