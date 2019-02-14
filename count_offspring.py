import sys
import phylonetwork
from generating_TC import count_feasible_pairs
import multiprocessing
import time


def treat_line(line):
    line = line.strip()
    net = phylonetwork.PhyloNetwork(eNewick=line)
    return count_feasible_pairs(net)


if __name__ == '__main__':
    try:
        filein = sys.argv[1]
    except IndexError:
        print("Usage: python count_offspring.py infile")
        exit()

    pool = multiprocessing.Pool()
    start = time.time()

    with open(filein, 'r') as f:
        lines = f.readlines()
        numlines = len(lines)
        outs = pool.imap_unordered(treat_line, lines)
        proclines = 0

        numnets = 0
        for out in outs:
            proclines += 1
            numnets += out
            now = time.time()
            print(f"Processed {proclines} of {len(lines)} lines. " +
                  f"Elapsed time: {now - start:.2f}s Expected time: {(now - start) * numlines / proclines:.2f}s " +
                  f"Network count: {numnets}. " +
                  f"Ratio: {numnets / proclines:.2f} Expected: {int(numnets * numlines / proclines)}")
        print(f"Total offpring: {numnets} networks")
