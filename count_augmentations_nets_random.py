import sys
import phylonetwork
from generating_TC import count_augmentations
import multiprocessing
import random
import math

try:
    filein = sys.argv[1]
except:
    filein = "nets_5.txt"

def treat_line(line):
    line = line.strip()
    net = phylonetwork.PhyloNetwork(eNewick=line)
    return count_augmentations(net)

if __name__ == '__main__':
    # pool = multiprocessing.Pool()

    with open(filein,'r') as f:
        lines = f.readlines()
        # outs = pool.imap_unordered(treat_line,lines)
        proclines = 0
        numnets = 0
        ratios = []
        while True:
            proclines += 1
            line = random.choice(lines)
            numnets += treat_line(line)
            newratio = numnets / proclines
            ratios.append(newratio)
            if proclines % 10 == 0:
                maxratio = max(ratios[-10:])
                minratio = min(ratios[-10:])
                if math.fabs((maxratio-minratio)/maxratio) < 0.001:
                    break
            print(f"In {proclines} of {len(lines)}. \tOut {numnets} \tRatio {newratio} \tEstimate: {newratio*len(lines)}")
        print(f"Last Estimate: {newratio*len(lines)}")
