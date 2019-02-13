import sys
import phylonetwork
from generating_TC import BTC_offspring_generator
import multiprocessing


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
        print("Usage: augmentate_nets.py infile outfile taxon")
        exit()

    pool = multiprocessing.Pool()

    with open(filein, 'r') as f, open(fileout, 'w') as g:
        lines = f.readlines()
        outs = pool.imap_unordered(treat_line, lines)
        numnets = 0
        for out in outs:
            for l in out:
                numnets += 1
                g.write(l+"\n")
            print(f"Currently {numnets}")
