"""
@author: gpwolfe
"""

import pandas as pd

def constrain(data1, data2)
data1 = pd.Series(data1)
data2 = pd.Series(data2)
bins = pd.IntervalIndex.from_tuples([(0, 0.35), (0.35, 0.7), (0.7, 1.1)],
                                    closed='right')
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