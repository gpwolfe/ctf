"""
@author: gpwolfe
"""
import os
import pandas as pd

def constrain(data1, data2, bins)
data1 = pd.Series(data1)
data2 = pd.Series(data2)

data1_comp = pd.cut(data1, bins=bins)
data2_comp = pd.cut(data2, bins=bins)
equal = data1_comp.eq(data2_comp)
reindex = equal.reindex(range(-2, len(equal)+2), fill_value=True)
window = reindex.rolling(5).sum().shift(-2).loc[0:len(equal)-1]
constraints = (window.values == 0) | (window.values == 5)

def load_files(directory):
    data = []
    for fn in sorted(os.listdir(directory)):
        if fn.endswith('.txt'):
            data.append(pd.read_table(os.path.join(directory, fn),
                                      header=None, index_col=0))
    return data

def find_constraints(directory):
    data = load_files(directory)
    
    bins = pd.IntervalIndex.from_tuples([(-.00001, 0.35), (0.35, 0.7),
                                         (0.7, 1.1)], closed='right')
    first = pd.cut(data[0].iloc[:,0], bins=bins)
    data1 = first.iloc[:]
    for struct in data[1:]:
        


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