#!/usr/bin/env python
from random import sample, expovariate
import argparse
import lineage
import pickle
import sys
def coalescence_time(i, n=1):
    return expovariate(n * (i * (i - 1)) / 2.0)

def mutation_time(theta, i):
    return expovariate(theta / 2 * i)

def coalesce(prev_time, lineages, n_mutations, t0=None, theta=None):
    """Select two lineages to merge after an appropriate period of time.
    
    ARGS:
        prev_time:  the amount of time that has elapsed prior to this
            coalescence event
        lineages:  the collection of lineages from which to select and coalesce

    RETURNS:
        prev_time plus the coalescent time for this coalescence
    """
    left, right = sample(lineages, 2)
    if t0 is not None:
        if prev_time > t0:
            # if we start over t0, use the smaller population
            time = prev_time + coalescence_time(len(lineages), n=1)
        else:
            coal_time = coalescence_time(len(lineages), n=2)
            if coal_time + prev_time < t0:
                # if we're still under t0, behave as normal
                time = coal_time + prev_time
            else:
                # halve the time beyond t0 to account for faster coalescence
                # with smaller population
                time = 0.5*(coal_time + prev_time - t0) + t0
    else:
        time = prev_time + coalescence_time(len(lineages))
    if theta is not None:
        mut_time = mutation_time(theta, len(lineages)) + prev_time
        # Mutation occurs first
        if mut_time < time:
            l = sample(lineages, 1)[0]
            l.mutate(n_mutations)
            return (mut_time, n_mutations + 1)
    new_lineage = lineage.Lineage(time, left, right)
    lineages.remove(left)
    lineages.remove(right)
    lineages.add(new_lineage)
    return (time, n_mutations)

def create_starting_population(i):
    """Create a starting population of i individuals."""
    lineages = set()
    for j in range(0, i):
        lineages.add(lineage.Lineage(0.0, None, None))
    return lineages

def simulate(n, t0=None, theta=None):
    """Run a simulation on a sample of size n"""
    lineages = create_starting_population(n)
    time = 0
    n_mutations = 0
    while len(lineages) > 1:
        (time, n_mutations) = coalesce(time, lineages, n_mutations,
                                       t0=t0, theta=theta)
    return lineages.pop(), n_mutations

def e_eta_i(theta, i, n):
    """Return expected value of eta_i for sample size n"""
    if i == n-1:
        delta = 1
    else:
        delta = 0
    return theta * ((1.0/i) + (1/(n-1)))/(1 + delta)

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
    parser.add_argument('-t', '--t0', default=None, help='time at which '
                        'growth from N to present 2N occurred')
    parser.add_argument('--theta', default=None, help='scaled mutation rate')
    parser.add_argument('--sequences', action='store_true', default=False)
    parser.add_argument('--frequency', action='store_true', default=False)
    args = parser.parse_args()

    if args.input:
        root = load(args.input)
    else:
        if args.t0 is not None:
            args.t0 = float(args.t0)
        if args.theta is not None:
            args.theta = float(args.theta)
        root, n_mutations = simulate(args.n, t0=args.t0, theta=args.theta)

    if args.output:
        save(root, args.output)

    branches = lineage.branch_length(root)

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
            print '\n'.join([str(x) for x in lineage.pairwise_times(root)])
    if args.sequences or args.frequency:
        sequences = lineage.build_sequences(root)
        converted = [lineage.convert_sequence(x,n_mutations) for x in sequences]
        if args.sequences:
            stringified = [[str(y) for y in x] for x in converted]
            print '\n'.join([''.join(x) for x in stringified])
        if args.frequency:
            frequency = lineage.build_site_frequency(converted)
            for i in range(1, max(frequency.keys() + [-1, args.n/2]) + 1):
                if i in frequency:
                    print "{}\t{:f}\t{:f}".format(i,
                                       float(frequency[i])/args.n,
                                       e_eta_i(args.theta, i, args.n))
                else:
                    print "{}\t{:f}\t{:f}".format(i, 0,
                                       e_eta_i(args.theta, i, args.n))

if __name__ == '__main__':
    main()
