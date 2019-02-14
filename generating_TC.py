#!/usr/bin/env python
# coding: utf-8

import phylonetwork
import networkx
from itertools import product, permutations, combinations, combinations_with_replacement
from networkx.drawing.nx_agraph import graphviz_layout
import random


def split(net, u, hybrid=False):
    """
    Splits the node u by introducing an elementary node. If hybrid is True, then
    it is assumed that it will be a hybrid node and hence its name has '#' as first
    character; otherwise it relies on the _generate_new_id() method.
    """
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


def is_feasible_common(net, tis, yis):
    """
    Tests if the common conditions for a feasible pair are satisfied
    for the pair (tis,yis) inside net.
    """
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
            if True not in [ti == yk for ti in tis for yk in [yi, yj]]:
                return False
    return True


def is_feasible_T(net, tis, yis):
    """
    Tests if the conditions for a T-feasible pair are satisfied
    for the pair (tis,yis) inside net.
    """
    if len(tis) != 1:
        return False
    return is_feasible_common(net, tis, yis)


def is_feasible_H(net, tis, yis):
    """
    Tests if the conditions for a H-feasible pair are satisfied
    for the pair (tis,yis) inside net.
    """
    if not len(tis) in [1, 2]:
        return False
    for ti in tis:
        if ti in yis:
            return False
    return is_feasible_common(net, tis, yis)


def augmentation_T(net, ell, tis, yis, pre_test=True):
    """
    Finds the T-augmentation network of net using the pair (tis,yis) and adding
    the taxon ell (if possible; otherwise returns None). If pre_test is True,
    it test if (tis,yis) for a T-feasible pair; otherwise it first computes
    the network and then checks if it is a BTC network.
    """
    if pre_test and not is_feasible_T(net, tis, yis):
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
    """
    Finds the H-augmentation network of net using the pair (tis,yis) and adding
    the taxon ell (if possible; otherwise returns None). If pre_test is True,
    it test if (tis,yis) for a H-feasible pair; otherwise it first computes
    the network and then checks if it is a BTC network.
    """
    if pre_test and not is_feasible_H(net, tis, yis):
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
    us = [u0]+[netb._generate_new_id() for _ in range(1, rt)]
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


def draw(net):
    """
    Draws the network net using dot layout.
    """
    pos = graphviz_layout(net, prog='dot')
    networkx.draw_networkx(net, pos, labels=net._labels)


def is_tree_child(net):
    """
    Tests if net is tree-child
    """
    for u in net.nodes():
        if net.is_leaf(u):
            continue
        children = net.successors(u)
        if True not in [net.is_tree_node(child) for child in children]:
            return False
    return True


def feasible_pairs_generator(net):
    """
    Generator that yields all feasible pairs inside net.
    Each yielded value has the form [augmentation_function, S1, S2],
    where augmentation_function is the function that has to be applied to
    net taking as recovering data S1 and S2.
    """
    nodes = net.nodes()
    tree_nodes = set([u for u in nodes if net.is_tree_node(u)])
    hyb_nodes = set([u for u in nodes if net.is_hybrid_node(u)])
    num_hyb_nodes = len(hyb_nodes)
    leaves = net.leaves()
    num_leaves = len(leaves)
    for r in range(num_leaves - num_hyb_nodes + 1):  # num_tree_nodes + 1):
        for (wt, vts) in product(tree_nodes, permutations(tree_nodes, r)):
            if is_feasible_T(net, {wt}, vts):
                yield [augmentation_T, {wt}, vts]
    for r in range(num_leaves - num_hyb_nodes):
        for ts in combinations_with_replacement(tree_nodes, 2):
            tsset = set(ts)
            for yis in permutations(tree_nodes - tsset, r):
                if is_feasible_H(net, tsset, yis):
                    yield [augmentation_H, tsset, yis]


def feasible_pairs(net):
    return list(feasible_pairs_generator(net))


def count_feasible_pairs(net):
    return len(feasible_pairs(net))


def BTC_offspring_generator(net, ell):
    """
    Generator that yields all the offspring of the network net having ell as
    the new taxon.
    """
    for feasible in feasible_pairs_generator(net):
        f = feasible[0]
        args = feasible[1:]
        yield f(net, ell, *args)


def BTC_networks_generator(taxa):
    """
    Generator that yields all the BTC networks over taxa.
    """
    if len(taxa) == 1:
        yield phylonetwork.PhyloNetwork(eNewick=taxa[0] + ';')
        return
    ell = taxa[-1]
    parent_generator = BTC_networks_generator(taxa[:-1])
    for net in parent_generator:
        for augmented in BTC_offspring_generator(net, ell):
            yield augmented


def random_BTC_network(taxa):
    """
    Returns a random BTC network over taxa. It uses a recursive implementation, where
    at each stage the generated network is chosen uniformly between its offspring
    (which does not imply that the overall procedure is uniform).
    """
    if len(taxa) == 1:
        return phylonetwork.PhyloNetwork(eNewick=taxa[0] + ';')
    net_previous = random_BTC_network(taxa[:-1])
    feasibles = feasible_pairs(net_previous)
    chosen = random.choice(feasibles)
    f = chosen[0]
    args = chosen[1:]
    return f(net_previous, taxa[-1], *args)
