import sys
import phylonetwork
from generating_TC import count_augmentations
import multiprocessing

filein = sys.argv[1]

def treat_line(line):
    line = line.strip()
    net = phylonetwork.PhyloNetwork(eNewick=line)
    return count_augmentations(net)

if __name__ == '__main__':
    pool = multiprocessing.Pool()

    with open(filein,'r') as f:
        lines = f.readlines()
        outs = pool.imap_unordered(treat_line,lines)
        proclines = 0

        numnets = 0
        for out in outs:
            proclines += 1
            numnets += out
            print(f"In {proclines} of {len(lines)}. Out {numnets} Ratio {numnets / proclines}")
