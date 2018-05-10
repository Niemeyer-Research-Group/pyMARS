from collections import deque
from heapq import heappush, heappop
from itertools import count
import networkx as nx
from networkx.utils import generate_unique_node

#The following functions were taken from the networkx source code and modified for this project. 
 
def ss_dijkstra_path_length_modified(G, source, cutoff=None,
                                       weight='weight'):
    """Compute the greatest path length via multiplication between source and all other
    reachable nodes for a weighted graph with all weights <= 1.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       Starting node for path

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight.

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    length : dictionary
       Dictionary of shortest lengths keyed by target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> length=nx.ss_dijkstra_path_length_modified(G,0)
    >>> length[4]
    1
    >>> print(length)
    {0: 0, 1: 1, 2: 1, 3: 1, 4: 1}

    Notes
    -----
    Edge weight attributes must be numerical and <= 1.
    Distances are calculated as products of weighted edges traversed.
    Don't use a cutoff.

    See Also
    --------
    single_source_dijkstra()

    """
    if G.is_multigraph():
        get_weight = lambda u, v, data: min(
            eattr.get(weight, 1) for eattr in data.values())
    else:
        get_weight = lambda u, v, data: data.get(weight, 1)

    return mod_dijkstra(G, source, get_weight, cutoff=cutoff)



def mod_dijkstra(G, source, get_weight, pred=None, paths=None, cutoff=None,
              target=None):
    """Modified implementation of Dijkstra's algorithm. Multiples values instead of adds and returns a dictionary with nodes as keys and values containing the greatest path to that node.  Each edge weight must be <= 1 so that the further they are from the source, the less important they are.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       Starting node for path

    get_weight: function
        Function for getting edge weight

    pred: list, optional(default=None)
        List of predecessors of a node

    paths: dict, optional (default=None)
        Path from the source to a target node.

    target : node label, optional
       Ending node for path

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    distance,path : dictionaries
       Returns a tuple of two dictionaries keyed by node.
       The first dictionary stores distance from the source.
       The second stores the path from the source to that node.

    pred,distance : dictionaries
       Returns two dictionaries representing a list of predecessors
       of a node and the distance to each node.

    distance : dictionary
       Dictionary of greatest lengths keyed by target.
    """
    G_succ = G.succ if G.is_directed() else G.adj #Adjaceny list
    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {source: 0}
    c = count()
    fringe = []  # use heapq with (distance,label) tuples
    push(fringe, (0, next(c), source))
    while fringe:
        cont = 0
        (d, _, v) = pop(fringe)
        if v in dist and d < dist[v]:
            continue  # already searched this node.
        if v == source:
            d = 1
        dist[v] = d 
        if v == target:
            break
        for u, e in G_succ[v].items(): #For all adjancent edges, get weights and multiply them by current path taken. 
            cost = get_weight(v, u, e)
            if cost is None:
                continue
            vu_dist = dist[v] * get_weight(v, u, e)
            if cutoff is not None:
                if vu_dist < cutoff:
                    continue
            #if v in dist:
                #if vu_dist > dist[v]:
                    #raise ValueError('Contradictory paths found:',
                                    #'weights greater than one?')
            elif u not in seen or vu_dist > seen[u]: #If this path to u is greater than any other path we've seen, push it to the heap to be added to dist.
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if paths is not None:
                    paths[u] = paths[v] + [u]
                if pred is not None:
                    pred[u] = [v]
            elif vu_dist == seen[u]:
                if pred is not None:
                    pred[u].append(v)

    if paths is not None:
        return (dist, paths)
    if pred is not None:
        return (pred, dist)
    return dist

