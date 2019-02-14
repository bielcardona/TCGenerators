import sys
import phylonetwork
from generating_TC import BTC_offspring_generator
import multiprocessing
import time


def treat_line(line):
    line = line.strip()
    net = phylonetwork.PhyloNetwork(eNewick=line)
    new_nets = BTC_offspring_generator(net, newlabel)
    return [net.eNewick() for net in new_nets]


if __name__ == '__main__':
    try:
        filein = sys.argv[1]
        fileout = sys.argv[2]
        newlabel = sys.argv[3]
    except IndexError:
        print("Usage: python compute_offspring.py infile outfile taxon")
        exit()

    pool = multiprocessing.Pool()
    start = time.time()

    with open(filein, 'r') as f, open(fileout, 'w') as g:
        lines = f.readlines()
        numlines = len(lines)
        outs = pool.imap_unordered(treat_line, lines)
        proclines = 0
        numnets = 0
        for out in outs:
            proclines += 1
            numnets += len(out)
            now = time.time()
            print(f"Processed {proclines} of {len(lines)} lines. " +
                  f"Elapsed time: {now - start:.2f}s Expected time: {(now - start) * numlines / proclines:.2f}s " +
                  f"Network count: {numnets}. " +
                  f"Ratio: {numnets / proclines:.2f} Expected: {int(numnets * numlines / proclines)}")
            for l in out:
                g.write(l + "\n")
        print(f"Total offpring: {numnets} networks")
