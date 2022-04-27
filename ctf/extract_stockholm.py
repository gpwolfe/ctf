import os

fns = sorted([fn for fn in os.listdir() if fn.endswith('out.txt')])
count = 1
for fn in fns:
    with open(f'CTF{count}_out.txt', 'r') as fin:
        data = fin.read().splitlines(True)
        stockholm_energy = data[2:]

        just_stockholm = stockholm_energy[0].split()[0]

        # print(just_stockholm)

    with open(f'CTF{count}_out_stockholm.txt', 'w') as fout:

        fout.writelines(just_stockholm)
    count += 1
    