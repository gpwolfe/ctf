"""
Generate RNA folding constraints based on similar reactivities.

Given a file containing known nucleotide bindings for an RNA fragment,
module returns a file containing probable preserved bindings for a lengthened
RNA fragment, based on the SHAPE reactivity for both fragments.

This can be used to model a longer RNA structure by iteratively preserving
folding constraints from shorter RNA fragments if the SHAPE reactivity does
not change by more than a defined amount.

Makes use of the parse_vienna_to_pairs function from rna-tools
(https://rna-tools.readthedocs.io).

SHAPE reactivity files should be saved using the following file name format:
"CTF{ctf number}.shape.txt" -->  CTF4.shape.txt

Example
-------
>>> python3 shape_constraints.py ctf2 ctf3

Output file
-----------
CTF3_constraints.txt

@author: gpwolfe
"""
from argparse import ArgumentParser
import sys

import pandas as pd


def constrain_bin(data1, data2):
    """
    Find indices of nucleotides that remain in same reactivity bin.

    Low, medium and high reactivity are defined as separate bins. If a
    nucleotide is in the same bin in both SHAPE files, this index is included
    in the returned Series.
    """
    equal = data1.eq(data2)
    reindex = equal.reindex(range(equal.index[0]-1, len(equal)+1),
                            fill_value=equal.iloc[0])
    reindex = reindex.reindex(range(reindex.index[0], reindex.index[-1]+2),
                              fill_value=reindex.iloc[-1])
    window = reindex.rolling(3).sum().shift(-1).loc[1:len(equal)]
    const_bin = window.values == 3
    return pd.Series(const_bin)


def constrain(shape1, shape2):
    """Find indices of nts with similar SHAPE reactivity across files."""
    shape1.loc[:, '2'] = shape2.iloc[:, 0]

    constraints = (
        (shape1.iloc[:, 0] >= shape1.iloc[:, 1] - shape1.iloc[:, 1] * .2)
        & (shape1.iloc[:, 0] <= shape1.iloc[:, 1] + shape1.iloc[:, 1] * .2)
            )
    return constraints


def find_constraints(ctf1, ctf2):
    """
    Generate set of indices for nucleotides with similar reactivity.

    Parameters
    ----------
    ctf1 : string
        'CTF' + number corresponding to shorter RNA fragment.
    ctf2 : string
        'CTF' + number corresponding to longer RNA fragment.

    Returns
    -------
    Set
        Combined indices of nucleotides with preserved reactivity.

    """
    shape1 = pd.read_table(f'{ctf1.upper()}.shape.txt', header=None,
                           index_col=0)
    shape2 = pd.read_table(f'{ctf2.upper()}.shape.txt', header=None,
                           index_col=0)

    bins = pd.IntervalIndex.from_tuples([(-.00001, 0.4), (0.4, 0.75),
                                         (0.75, 1)], closed='right')
    data1 = pd.cut(shape1.iloc[:, 0], bins=bins).iloc[:]
    data2 = pd.cut(shape2.iloc[:, 0], bins=bins).iloc[:]
    const_bin = constrain_bin(data1, data2)
    const_bin_ix = set(const_bin[const_bin].index)
    constraints = constrain(shape1, shape2)
    constraints_index = set(constraints[constraints].index)

    return constraints_index | const_bin_ix


def extract_stockholm(ctf1, ctf2):
    """
    Generate constraint files.

    Constraints generated based on two SHAPE reactivity files and one
    structure file in Stockholm format.

    Parameters
    ----------
    ctf1 : string
        'CTF' + number corresponding to SHAPE file and Stockholm file of
        shorter RNA fragment.
    ctf2 : string
        'CTF' + number corresponding to SHAPE file of longer RNA fragment.

    Returns
    -------
    None.

    """
    with open(f'{ctf1}_out_stockholm.txt', 'r') as fin:
        data = fin.read()

    parsed = pd.Series(parse_vienna_to_pairs(data)[0])
    constraints_ix = find_constraints(ctf1, ctf2)

    const_pairs = []
    for pair in parsed:
        if (pair[0] in constraints_ix) and (pair[1] in constraints_ix):
            const_pairs.append(pair)

    towrite = ""
    for pair in const_pairs:
        towrite += f'{pair[0]} {pair[1]}\n'
    with open(f'{ctf2.upper()}_constraints.txt', 'w') as fout:
        the_header = "DS:\n-1\nSS:\n-1\nMod:\n-1\nPairs:\n"
        the_footer = "-1 -1\nFMN:\n-1\nForbids:\n-1 -1"
        fout.writelines(the_header + towrite + the_footer)


class ExceptionOpenPairsProblem(Exception):
    """Exception passed to parse_vienna_to_pairs."""

    pass


def parse_vienna_to_pairs(ss, remove_gaps_in_ss=False):
    """Parse Vienna (dot-bracket notation) to get pairs.

    Args:

       ss (str): secondary stucture in Vienna (dot-bracket notation) notation
       remove_gaps_in_ss (bool): remove - from ss or not, design for DCA (tpp case
                                 ``ss = "(((((((((.((((.(((.....))))))......------)....."``
                                 works with pk of the first level, ``[[]]``

    Returns
    -------
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
        raise ExceptionOpenPairsProblem(
            'Too many open pairs (()) in structure')
    if stack_pk:
        raise ExceptionOpenPairsProblem(
            'Too many open pairs [[]] in structure')

    pairs.sort()
    pairs_pk.sort()
    return(pairs, pairs_pk)


def cmdline_exec(argv):
    """Generate binding constraints for RNA folding."""
    parser = ArgumentParser(
        description='Generate binding constraints for RNA folding')
    parser.add_argument('ctf1', nargs='?',
                        help='i.e.: "CTF1".')
    parser.add_argument('ctf2', help='i.e.: CTF2')

    args = parser.parse_args(argv)
    ctf1 = args.ctf1
    ctf2 = args.ctf2
    extract_stockholm(ctf1, ctf2)


if __name__ == '__main__':
    cmdline_exec(sys.argv[1:])
