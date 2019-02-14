import sys
import phylonetwork
from generating_TC import count_feasible_pairs
import multiprocessing


def treat_line(line):
    line = line.strip()
    net = phylonetwork.PhyloNetwork(eNewick=line)
    return count_feasible_pairs(net)


if __name__ == '__main__':
    try:
        filein = sys.argv[1]
    except IndexError:
        print("Usage: count_augmentation_nets infile")
        exit()

    pool = multiprocessing.Pool()

    with open(filein,'r') as f:
        lines = f.readlines()
        outs = pool.imap_unordered(treat_line,lines)
        proclines = 0

        numnets = 0
        for out in outs:
            proclines += 1
            numnets += out
            print(f"Processed {proclines} of {len(lines)} lines. Network count: {numnets}. Ratio: {numnets / proclines:.2f}")
