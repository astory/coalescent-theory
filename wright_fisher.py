#!/usr/bin/env python
from collections import Counter, namedtuple
from random import random, sample
import argparse

Gene = namedtuple('Gene', ['allele'])

def pick_new_generation(old_generation):
    """Make a new generation choosing randomly from the old, with replacement
    
    Note that this function does not copy the existing objects; the new list
    contains the *same* objects as before.
    """
    new_generation = []
    for _ in old_generation:
        new_generation.append(sample(old_generation, 1)[0])
    return new_generation

def create_starting_population(n, p):
    """Create a random diallelic population of size n

    ARGS:
        n: number of individuals in the population
        p: probability of allele 1
    """
    a1 = Gene(1)
    a2 = Gene(2)
    population = []
    for i in range(n):
        if random() < p:
            population.append(a1)
        else:
            population.append(a2)
    return population

def fix_starting_population(n, p):
    """Create a diallelic population of size n with fixed ratios

    ARGS:
        n: number of individuals in the population
        p: allele fraction of first allele
    """
    a1 = Gene(1)
    a2 = Gene(2)
    population = []
    for i in range(n):
        if i <= p * n:
            population.append(a1)
        else:
            population.append(a2)
    return population

def heterozygosity(population):
    """Calc the heterozygosity of a population, H=2p(1-p) with two alleles
    
    If more than two alleles are present, calculate heterozygosity assuming that
    all alleles beyond the first are the same.
    """
    c = Counter(population)
    (_, n) = c.most_common(1)[0]
    p = float(n) / len(population)
    return 2 * p * (1 - p)

def wright_fisher(n=10, gens=20, p=0.3, fix=False):
    """Run a wright_fisher model, return list of heterozygosities at gens"""
    if fix:
        pop = fix_starting_population(n, p)
    else:
        pop = create_starting_population(n, p)
    hets = [heterozygosity(pop)]
    for i in range(gens):
        pop = pick_new_generation(pop)
        hets.append(heterozygosity(pop))
    return hets

def main():
    parser = argparse.ArgumentParser(
                        description='Run a Wright-Fisher simulation')
    parser.add_argument('-f', '--fix', action='store_true', default=False)
    parser.add_argument('-n', '--population-size', default=10, type=int,
                        dest='n')
    parser.add_argument('-g', '--generation-number', default=20, type=int,
                        dest='gens')
    parser.add_argument('-p', '--allele-probability', default=0.3, type=float,
                        dest='p')
    parser.add_argument('-t', '--trials', default=1, type=int,
                        dest='trials')
    args = parser.parse_args()

    # Run a number of trials, and arrange them in columns for graphing
    trials = zip(range(args.gens + 1),
                 *[wright_fisher(n=args.n, gens=args.gens,
                                 p=args.p, fix=args.fix)
                   for _ in range(args.trials)])
    # Pretty print the trial results so they're easy to graph
    formatted = "\n".join(
        ["\t".join(
            [ '{:f}'.format(h) for h in trial])
        for trial in trials])
    print formatted

if __name__ == '__main__':
    main()
