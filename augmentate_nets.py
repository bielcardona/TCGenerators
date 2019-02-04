import sys
import phylonetwork
from generating_TC import find_augmentations
import multiprocessing

filein = sys.argv[1]
fileout = sys.argv[2]
newlabel = sys.argv[3]

def treat_line(line):
    line = line.strip()
    net = phylonetwork.PhyloNetwork(eNewick=line)
    new_nets = find_augmentations(net, newlabel)
    return [net.eNewick() for net in new_nets]

if __name__ == '__main__':
    pool = multiprocessing.Pool()

    with open(filein,'r') as f, open(fileout,'w') as g:
        lines = f.readlines()
        outs = pool.imap_unordered(treat_line,lines)
        numnets = 0
        for out in outs:
            for l in out:
                numnets += 1
                g.write(l+"\n")
            print(f"Currently {numnets}")
