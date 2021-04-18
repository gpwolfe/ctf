"""
@author: gpwolfe
"""
import os
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


def find_constraints(directory):
    constraints_index = []
    shape_data = load_shape(directory)
    bins = pd.IntervalIndex.from_tuples([(-.00001, 0.35), (0.35, 0.7),
                                         (0.7, 1.1)], closed='right')
    data1 = pd.cut(shape_data[0].iloc[:, 0], bins=bins).iloc[:]

    for struct in shape_data[1:]:
        data2 = pd.cut(struct.iloc[:, 0], bins=bins).iloc[:]
        
        constraints = pd.Series(constrain(data1, data2),
                                index=struct.index)
        constraints_index.append(set(constraints[constraints].index))
        data1 = data2.loc[:]
    return constraints_index


def extract_stockholm(directory):
    stock_fns = sorted([fn for fn in os.listdir() if fn.endswith('out.txt')])
    count1 = 1
    parsed = []
    for fn in stock_fns:
        with open(f'CTF{count1}_out.txt', 'r') as fin:
            data = fin.read().splitlines(True)
            stockholm_energy = data[2:]

            just_stockholm = stockholm_energy[0].split()[0]

            parsed.append(pd.Series(parse_vienna_to_pairs(just_stockholm)[0]))
            count1 += 1

    constraints_ix = find_constraints(directory)
    shape_stock = list(zip(constraints_ix, parsed[1:]))

    const_sets = {}
    count2 = 2
    for line in shape_stock:
        const_pairs = []
        for pair in line[1]:
            if (pair[0] in line[0]) and (pair[1] in line[0]):
                const_pairs.append(pair)
        const_sets[count2] = const_pairs
        count2 += 1

    for ix, key in enumerate(const_sets):
        with open(f'CTF{key}_constraints.txt', 'w') as fout:
            for pair in const_sets[key]:
                towrite = f'{pair[0]}, {pair[1]}\n'
                fout.writelines(towrite)
            




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