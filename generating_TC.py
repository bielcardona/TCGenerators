#!/usr/bin/env python
# coding: utf-8

import phylonetwork
import networkx
from itertools import product, permutations, combinations, combinations_with_replacement
from networkx.drawing.nx_agraph import graphviz_layout


def split(net, u, hybrid=False):
    if hybrid:
        h = len([u for u in net.nodes() if net.is_hybrid_node(u)])
        utilde = "#" + str(h+1)
    else:
        utilde = net._generate_new_id()
    for e in net.in_edges(u):
        net.add_edge(e[0], utilde)
        net.remove_edge(e[0], u)
    net.add_edge(utilde, u)
    return utilde


def is_correct_common(net, tis, yis):
    for ti in tis:
        if not net.is_tree_node(ti):
            return False
    for yi in yis:
        if not net.is_tree_node(yi):
            return False
        parents = net.predecessors(yi)
        if not parents:
            continue
        if yi in tis:
            continue
        parent = parents[0]
        if net.is_hybrid_node(parent):
            return False
        for child in net.successors(parent):
            if net.is_hybrid_node(child):
                return False
    for yi in yis:
        for ti in tis:
            if ti in networkx.descendants(net, yi):
                return False
    for (yi, yj) in combinations(yis, 2):
        if yi == yj:
            return False
        if net.predecessors(yi) == net.predecessors(yj):
            if not True in [ti == yk for ti in tis for yk in [yi, yj]]:
                return False
    return True


def is_correct_T(net, tis, yis):
    if len(tis) != 1:
        return False
    return is_correct_common(net, tis, yis)


def is_correct_H(net, tis, yis):
    if not len(tis) in [1, 2]:
        return False
    for ti in tis:
        if ti in yis:
            return False
    return is_correct_common(net, tis, yis)


def augmentation_T(net, ell, tis, yis, pre_test=True):
    if pre_test and not is_correct_T(net, tis, yis):
        return None
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
        v = split(netb, yis[i], hybrid=True)
        netb.add_edge(u, v)
    if not pre_test:
        if networkx.is_directed_acyclic_graph(netb) and is_tree_child(netb):
            return netb
        else:
            return None
    return netb


def augmentation_H(net, ell, tis, yis, pre_test=True):
    if pre_test and not is_correct_H(net, tis, yis):
        return None
    if len(tis) == 1:
        t1 = list(tis)[0]
        t2 = t1
    else:
        t1 = list(tis)[0]
        t2 = list(tis)[1]
    netb = net.copy()
    netb.cache = {}
    rt = len(yis) + 2
    h = len([u for u in netb.nodes() if net.is_hybrid_node(u)])
    u0 = "#" + str(h+1)
    us = [u0]+[netb._generate_new_id() for _ in range(1,rt)]
    netb.add_path(us)
    netb._labels[us[-1]] = ell
    w1 = split(netb, t1)
    netb.add_edge(w1, us[0])
    w2 = split(netb, t2)
    netb.add_edge(w2, us[0])
    for i in range(1, rt - 1):
        u = us[i]
        v = split(netb, yis[i - 1], hybrid=True)
        netb.add_edge(u, v)
    if not pre_test:
        if networkx.is_directed_acyclic_graph(netb) and is_tree_child(netb):
            return netb
        else:
            return None
    return netb


def BTC_networks(taxa):
    if len(taxa) == 1:
        yield phylonetwork.PhyloNetwork(eNewick=taxa[0] + ';')
        return
    ell = taxa[-1]
    parent_generator = BTC_networks(taxa[:-1])
    for net in parent_generator:
        nodes = net.nodes()
        tree_nodes = set([u for u in nodes if net.is_tree_node(u)])
        hyb_nodes = set([u for u in nodes if net.is_hybrid_node(u)])
        num_hyb_nodes = len(hyb_nodes)
        leaves = net.leaves()
        num_leaves = len(leaves)
        for r in range(num_leaves - num_hyb_nodes + 1):  # num_tree_nodes + 1):
            for (wt, vts) in product(tree_nodes, permutations(tree_nodes, r)):
                netb = augmentation_T(net, ell, {wt}, vts)
                if netb:
                    yield netb
        for r in range(num_leaves - num_hyb_nodes):
            for ts in combinations_with_replacement(tree_nodes, 2):
                tsset = set(ts)
                for yis in permutations(tree_nodes - tsset, r):
                    netb = augmentation_H(net, ell, tsset, yis)
                    if netb:
                        yield netb


def draw(net):
    pos = graphviz_layout(net, prog='dot')
    networkx.draw_networkx(net, pos, labels=net._labels)


def is_tree_child(net):
    for u in net.nodes():
        if net.is_leaf(u):
            continue
        children = net.successors(u)
        if True not in [net.is_tree_node(child) for child in children]:
            return False
    return True

# deprecated
def find_augmentations(net, ell, pre_test=True):
    nets = []
    nodes = net.nodes()
    tree_nodes = set([u for u in nodes if net.is_tree_node(u)])
    num_tree_nodes = len(tree_nodes)
    hyb_nodes = set([u for u in nodes if net.is_hybrid_node(u)])
    num_hyb_nodes = len(hyb_nodes)
    leaves = net.leaves()
    num_leaves = len(leaves)
    maxr_t = 0
    maxr_h = 0
    for r in range(num_leaves - num_hyb_nodes + 1):  # num_tree_nodes + 1):
        for (wt, vts) in product(tree_nodes, permutations(tree_nodes, r)):
            netb = augmentation_T(net, ell, {wt}, vts, pre_test)
            if netb:
                nets.append(netb)
                maxr_t = r
    for r in range(num_leaves - num_hyb_nodes):
        for ts in combinations_with_replacement(tree_nodes, 2):
            tsset = set(ts)
            for yis in permutations(tree_nodes - tsset, r):
                netb = augmentation_H(net, ell, tsset, yis, pre_test)
                if netb:
                    # print(tsset, yis)
                    maxr_h = r
                    nets.append(netb)
    print(f"n={num_leaves}, t={num_tree_nodes}, h={num_hyb_nodes}, rt={maxr_t}, rh={maxr_h}, #nets={len(nets)}")
    return nets


#deprecated
def find_networks(taxa, pre_test=True):
    if len(taxa) == 1:
        return [phylonetwork.PhyloNetwork(eNewick=taxa[0] + ';')]
    nets = []
    nets_previous = find_networks(taxa[:-1])
    num_nets = len(nets_previous)
    for i in range(num_nets):
        print(f"Taxa: {taxa}, {i} of {num_nets}")
        net = nets_previous[i]
        nets.extend(find_augmentations(net, taxa[-1], pre_test))
        print(len(nets))
    return nets
