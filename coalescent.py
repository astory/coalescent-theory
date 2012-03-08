#!/usr/bin/env python
from random import sample, expovariate
import argparse
import pickle

# This has to be a class rather than a named tuple because we want objects with
# identical contents (i.e., the starting lineages) to be distinct.
class Lineage():
    def __init__(self, time, left_child, right_child):
        # time since the start of the simulation when this lineage coalesced
        self.time = time
        self.left = left_child
        self.right = right_child

    def __str__(self):
        return "Lineage coalescing at " + str(self.time)

def coalescence_time(i):
    return expovariate((i * (i - 1)) / 2.0)

def coalesce(prev_time, lineages):
    """Select two lineages to merge after an appropriate period of time.
    
    ARGS:
        prev_time:  the amount of time that has elapsed prior to this
            coalescence event
        lineages:  the collection of lineages from which to select and coalesce

    RETURNS:
        prev_time plus the coalescent time for this coalescence
    """
    left, right = sample(lineages, 2)
    time = prev_time + coalescence_time(len(lineages))
    new_lineage = Lineage(time, left, right)
    lineages.remove(left)
    lineages.remove(right)
    lineages.add(new_lineage)
    return time

def create_starting_population(i):
    """Create a starting population of i individuals."""
    lineages = set()
    for j in range(0, i):
        lineages.add(Lineage(0.0, None, None))
    return lineages

def simulate(n):
    """Run a simulation on a sample of size n"""
    lineages = create_starting_population(n)
    time = 0
    while len(lineages) > 1:
        time = coalesce(time, lineages)
    return lineages.pop()

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

def save(tree, filehandle):
    pickle.dump(tree, filehandle)

def load(filehandle):
    return pickle.load(filehandle)

def main():
    parser = argparse.ArgumentParser( description='Run a coalescent simulation')
    parser.add_argument('-n', '--individuals', default=10, type=int,
                        dest='n')
    parser.add_argument('-q', '--quiet', action='store_true', default=False)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=None)
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=None)
    parser.add_argument('--T_MRCA', action='store_true', default=False)
    parser.add_argument('--branch_length', action='store_true', default=False)
    parser.add_argument('--tij', action='store_true', default=False)
    args = parser.parse_args()

    if args.input:
        root = load(args.input)
    else:
        root = simulate(args.n)

    if args.output:
        save(root, args.output)

    branches = branch_length(root)

    mrca = ''
    bl = ''
    if args.T_MRCA:
        mrca = str(root.time)
        if not args.quiet:
            mrca = 'T_MRCA ' + mrca
    if args.branch_length:
        bl = str(branches)
        if not args.quiet:
            bl = 'Branch Length: ' + bl
    else:
        if mrca or bl:
            print '\t'.join([mrca, bl])
        if args.tij:
            print '\n'.join([str(x) for x in pairwise_times(root)])

if __name__ == '__main__':
    main()
