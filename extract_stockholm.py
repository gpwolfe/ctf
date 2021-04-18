with open('CTF1_out.txt', 'r') as fin:
    data = fin.read().splitlines(True)
    stockholm_energy = data[2:]

    just_stockholm = stockholm_energy[0].split()[0]

    #print(just_stockholm)

with open('CTF1_out_stockholm.txt', 'w') as fout:

    fout.writelines(just_stockholm)