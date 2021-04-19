"""
@author: gpwolfe
"""
from argparser import ArgumentParser
import re
import sys

CT_HEADER_RE = re.compile(
    r'\s+\d+\s+ENERGY\s=\s(?P<energy>-?\d+\.\d+)\s+(?P<ctf>CTF\d+)')
CT_LINE_RE = re.compile(
    r'\s+\d+\s\w\s+\d+\s+\d+\s+(?P<start>\d*[1-9]\d*)\s+(?P<end>\d*[1-9]\d*)')


def parse_ct(fn):
    with open(fn) as fin:
        for line in fin:
            ln = line
            if CT_HEADER_RE.match(ln):
                head = CT_HEADER_RE.match(ln)
                energy = head.group('energy')
                ctf_num = head.group('ctf')
            else:
                with open(f'{ctf_num}_{energy}_cs.txt', 'a') as fout:
                    if CT_LINE_RE.match(ln):
                        ct_line = CT_LINE_RE.match(ln)
                        start = ct_line.group('start')
                        end = ct_line.group('end')
                        fout.write(f'{start} {end}\n')


def cmdline_exec(argv):

    parser = ArgumentParser(description='Parse multi-ctf data files.')
    parser.add_argument('filename', help='Path to file containing ctf data')
    args = parser.parse_args(argv)
    filename = args.filename
    parse_ct(filename)


if __name__ == '__main__':
    cmdline_exec(sys.argv[1:])