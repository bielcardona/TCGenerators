#!/usr/bin/env python
# coding: utf-8

import phylonetwork
import networkx
from itertools import product, permutations, combinations, combinations_with_replacement
from networkx.drawing.nx_agraph import graphviz_layout


def split(net, u):
    utilde = net._generate_new_id()
    for e in net.in_edges(u):
        net.add_edge(e[0], utilde)
        net.remove_edge(e[0], u)
    net.add_edge(utilde, u)
    return utilde


def is_correct_common(net, tis, yis):
    for ti in tis:
        if not net.is_tree_node(ti): return False
    for yi in yis:
        if not net.is_tree_node(yi): return False
        parents = net.predecessors(yi)
        if not parents: continue
        parent = parents[0]
        if net.is_hybrid_node(parent): return False
        for child in net.successors(parent):
            if net.is_hybrid_node(child): return False
        for ti in tis:
            if ti in networkx.descendants(net,yi): return False
    for (yi, yj) in combinations(yis, 2):
        if yi == yj: return False
        if net.predecessors(yi) == net.predecessors(yj): return False
    return True


def is_correct_T(net, tis, yis):
    if len(tis) != 1: return False
    return is_correct_common(net, tis, yis)


def is_correct_H(net, tis, yis):
    if not len(tis) in [1, 2]: return False
    for ti in tis:
        if ti in yis: return False
    return is_correct_common(net, tis, yis)


def augmentation_T(net, ell, tis, yis):
    if not is_correct_T(net, tis, yis): return None
    t1 = list(tis)[0]
    netb = net.copy()
    netb.cache = {}
    r = len(yis) + 1
    us = [netb._generate_new_id() for _ in range(r)]
    netb.add_path(us)
    netb._labels[us[-1]] = ell
    w = split(netb, t1)
    netb.add_edge(w, us[0])
    for i in range(r - 1):
        u = us[i]
        v = split(netb, yis[i])
        netb.add_edge(u, v)
    return netb


def augmentation_H(net, ell, tis, yis):
    if not is_correct_H(net, tis, yis): return None
    if len(tis) == 1:
        t1 = list(tis)[0]
        t2 = t1
    else:
        t1 = list(tis)[0]
        t2 = list(tis)[1]
    netb = net.copy()
    netb.cache = {}
    rt = len(yis) + 2
    us = [netb._generate_new_id() for _ in range(rt)]
    netb.add_path(us)
    netb._labels[us[-1]] = ell
    w1 = split(netb, t1)
    netb.add_edge(w1, us[0])
    w2 = split(netb, t2)
    netb.add_edge(w2, us[0])
    for i in range(1, rt - 1):
        u = us[i]
        v = split(netb, yis[i - 1])
        netb.add_edge(u, v)
    return netb


def find_augmentations(net, ell):
    nets = []
    nodes = net.nodes()
    tree_nodes = set([u for u in nodes if net.is_tree_node(u)])
    num_tree_nodes = len(tree_nodes)
    for r in range(num_tree_nodes + 1):
        for (wt, vts) in product(tree_nodes, permutations(tree_nodes, r)):
            netb = augmentation_T(net, ell, set([wt]), vts)
            if netb:
                print(wt, vts)
                nets.append(netb)
        for ts in combinations_with_replacement(tree_nodes, 2):
            tsset = set(ts)
            for yis in permutations(tree_nodes - tsset, r):
                netb = augmentation_H(net, ell, tsset, yis)
                if netb:
                    print(tsset, yis)
                    nets.append(netb)
    return nets


def find_networks(taxa):
    if len(taxa) == 1:
        return [phylonetwork.PhyloNetwork(eNewick=taxa[0] + ';')]
    nets = []
    for net in find_networks(taxa[:-1]):
        nets.extend(find_augmentations(net, taxa[-1]))
        print(len(nets))
    return nets


def draw(net):
    pos = graphviz_layout(net, prog='dot')
    networkx.draw_networkx(net, pos, labels=net._labels)

def is_tree_child(net):
    for u in net.nodes():
        if net.is_leaf(u): continue
        children = net.successors(u)
        if True not in [net.is_tree_node(child) for child in children]:
            return False
    return True
