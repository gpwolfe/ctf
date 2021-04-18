"""
@author: gpwolfe
"""
from argparse import ArgumentParser
import os
import sys

import pandas as pd


def constrain(data1, data2):
    equal = data1.eq(data2)
    reindex = equal.reindex(range(-2, len(equal)+2), fill_value=True)
    window = reindex.rolling(5).sum().shift(-2).loc[0:len(equal)-1]
    constraints = (window.values == 0) | (window.values == 5)
    return constraints

def load_shape(directory):
    shape_data = []
    for fn in sorted(os.listdir(directory)):
        if fn.endswith('shape.txt'):
            shape_data.append(pd.read_table(os.path.join(directory, fn),
                                            header=None, index_col=0))
    return shape_data


def find_constraints(ctf1, ctf2):

    shape1 = pd.read_table(f'{ctf1.upper()}.shape.txt', header=None,
                           index_col=0)
    shape2 = pd.read_table(f'{ctf2.upper()}.shape.txt', header=None,
                           index_col=0)

    bins = pd.IntervalIndex.from_tuples([(-.00001, 0.35), (0.35, 0.7),
                                         (0.7, 1.1)], closed='right')
    data1 = pd.cut(shape1.iloc[:, 0], bins=bins).iloc[:]
    data2 = pd.cut(shape2.iloc[:, 0], bins=bins).iloc[:]

    constraints = pd.Series(constrain(data1, data2),
                            index=data2.index)
    constraints_index = set(constraints[constraints].index)

    return constraints_index


def extract_stockholm(ctf1, ctf2):

    with open(f'{ctf1}_out.txt', 'r') as fin:
        data = fin.read().splitlines(True)
        stockholm_energy = data[2:]
        just_stockholm = stockholm_energy[0].split()[0]
        parsed = (pd.Series(parse_vienna_to_pairs(just_stockholm)[0]))

    constraints_ix = find_constraints(ctf1, ctf2)

    const_pairs = []
    for pair in parsed:
        if (pair[0] in constraints_ix) and (pair[1] in constraints_ix):
            const_pairs.append(pair)

    towrite = ""
    for pair in const_pairs:
        towrite += f'{pair[0]}, {pair[1]}\n'
    with open(f'{ctf2.upper()}_constraints.txt', 'w') as fout:
        the_header = "DS:\n-1\nSS:\n-1\nMod:\n-1\nPairs:\n"
        the_footer = "-1 -1\nFMN:\n-1\nForbids:\n-1 -1"
        fout.writelines(the_header + towrite + the_footer)


class ExceptionOpenPairsProblem(Exception):
    pass


def parse_vienna_to_pairs(ss, remove_gaps_in_ss=False):
    """Parse Vienna (dot-bracket notation) to get pairs.

    Args:

       ss (str): secondary stucture in Vienna (dot-bracket notation) notation
       remove_gaps_in_ss (bool): remove - from ss or not, design for DCA (tpp case
                                 ``ss = "(((((((((.((((.(((.....))))))......------)....."``
                                 works with pk of the first level, ``[[]]``

    Returns:

        list of two lists: (pairs, pairs_pk)

    Examples::

        >>> parse_vienna_to_pairs('((..))')
        ([[1, 6], [2, 5]], [])

        >>> parse_vienna_to_pairs('(([[))]]')
        ([[1, 6], [2, 5]], [[3, 8], [4, 7]])

        >>> parse_vienna_to_pairs('((--))')
        ([[1, 6], [2, 5]], [])

        >>> parse_vienna_to_pairs('((--))', remove_gaps_in_ss=True)
        ([[1, 4], [2, 3]], [])

        >>> parse_vienna_to_pairs('((((......')
        Traceback (most recent call last):
          File "/usr/lib/python2.7/doctest.py", line 1315, in __run
            compileflags, 1) in test.globs
          File "<doctest __main__.parse_vienna_to_pairs[4]>", line 1, in <module>
            parse_vienna_to_pairs('((((......')
          File "./SecondaryStructure.py", line 106, in parse_vienna_to_pairs
            raise ExceptionOpenPairsProblem('Too many open pairs (()) in structure')
        ExceptionOpenPairsProblem: Too many open pairs (()) in structure

    """
    if remove_gaps_in_ss:
        ss = ss.replace('-', '')
    stack = []
    pairs = []
    pairs_pk = []
    stack_pk = []
    for c, s in enumerate(ss):
        if s == '(':
            stack.append(c + 1)
        if s == ')':
            pairs.append([stack.pop(), c + 1])
        if s == '[':
            stack_pk.append(c + 1)
        if s == ']':
            pairs_pk.append([stack_pk.pop(), c + 1])

    if stack:
        raise ExceptionOpenPairsProblem('Too many open pairs (()) in structure')
    if stack_pk:
        raise ExceptionOpenPairsProblem('Too many open pairs [[]] in structure')

    pairs.sort()
    pairs_pk.sort()
    return(pairs, pairs_pk)


def cmdline_exec(argv):
    """Run electrostatics histrogram plotting from command line."""
    parser = ArgumentParser(
        description='')
    parser.add_argument('ctf1', nargs='?',
                        help='I.E.: "CTF1".')
    parser.add_argument('ctf2', help='Path to new file.')

    args = parser.parse_args(argv)
    ctf1 = args.ctf1
    ctf2 = args.ctf2
    extract_stockholm(ctf1, ctf2)

if __name__ == '__main__':
    cmdline_exec(sys.argv[1:])
