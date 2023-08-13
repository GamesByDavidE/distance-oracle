import collections
import ctypes
import math
import itertools
import random

import networkx


_distance_oracle = ctypes.CDLL("distance_oracle.dylib")

_NodeIndex = ctypes.c_uint32
_ArcIndex = ctypes.c_uint32
_Coordinate = ctypes.c_double
_Weight = ctypes.c_double
_DistanceOracleHandle = ctypes.c_void_p

_NewDistanceOracle = _distance_oracle.NewDistanceOracle
_NewDistanceOracle.restype = None
_NewDistanceOracle.argtypes = [
    ctypes.POINTER(_DistanceOracleHandle),  # graph_handle
    _NodeIndex,  # node_count
    ctypes.POINTER(_ArcIndex),  # first_outgoing_arcs
    ctypes.POINTER(_NodeIndex),  # arc_heads
]

_PlaceNode = _distance_oracle.PlaceNode
_PlaceNode.restype = None
_PlaceNode.argtypes = [
    _DistanceOracleHandle,  # graph_handle
    _NodeIndex,  # node
    _Coordinate,  # x
    _Coordinate,  # y
    _Coordinate,  # z
]

_LightenArc = _distance_oracle.LightenArc
_LightenArc.restype = None
_LightenArc.argtypes = [
    _DistanceOracleHandle,  # graph_handle
    _ArcIndex,  # arc
    _Weight,  # new_weight
]

_FindShortestPath = _distance_oracle.FindShortestPath
_FindShortestPath.restype = ctypes.c_bool
_FindShortestPath.argtypes = [
    _DistanceOracleHandle,  # graph_handle
    _NodeIndex,  # s
    _NodeIndex,  # t
    ctypes.POINTER(_NodeIndex),  # path_node_count
    ctypes.POINTER(_NodeIndex),  # path_nodes
]

_DeleteDistanceOracle = _distance_oracle.DeleteDistanceOracle
_DeleteDistanceOracle.restype = None
_DeleteDistanceOracle.argtypes = [
    ctypes.POINTER(_DistanceOracleHandle),  # graph_handle
]


class DistanceOracle:
    def __init__(self, nx_g):
        self._handle = _DistanceOracleHandle()
        self._nx_nodes = list(nx_g.nodes())
        self._node_indexes = {u: i for (i, u) in enumerate(self._nx_nodes)}
        self._arc_indexes = {}
        first_outgoing_arcs = (_ArcIndex * (nx_g.number_of_nodes() + 1))()
        arc_heads = (_NodeIndex * (2 * nx_g.number_of_edges()))()
        j = 0
        for i, u in enumerate(self._nx_nodes):
            for v in nx_g.neighbors(u):
                self._arc_indexes[(u, v)] = j
                arc_heads[j] = self._node_indexes[v]
                j += 1
            first_outgoing_arcs[i + 1] = j
        _NewDistanceOracle(
            ctypes.byref(self._handle),
            nx_g.number_of_nodes(),
            first_outgoing_arcs,
            arc_heads,
        )
        for i, u in enumerate(self._nx_nodes):
            _PlaceNode(self._handle, i, u.x, u.y, u.z)
        for u, v in nx_g.edges():
            self.lighten_edge(u, v, nx_g.get_edge_data(u, v)["weight"])

    def lighten_edge(self, u, v, new_weight):
        _LightenArc(self._handle, self._arc_indexes[(u, v)], new_weight)
        _LightenArc(self._handle, self._arc_indexes[(v, u)], new_weight)

    def find_shortest_path(self, s, t):
        s = self._node_indexes[s]
        t = self._node_indexes[t]
        path_node_count = _NodeIndex()
        path_nodes = (_NodeIndex * len(self._nx_nodes))()
        if not _FindShortestPath(
            self._handle,
            s,
            t,
            ctypes.byref(path_node_count),
            path_nodes,
        ):
            return None
        return [self._nx_nodes[u] for u in path_nodes[: path_node_count.value]]

    def __del__(self):
        _DeleteDistanceOracle(self._handle)


_Point = collections.namedtuple("_Point", ["x", "y", "z"])


def _distance(a, b):
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _total_distance(path):
    return (
        None
        if path is None
        else sum(_distance(a, b) for (a, b) in itertools.pairwise(path))
    )


def _test():
    nx_g = networkx.Graph()
    n = random.randrange(2, 100)
    while nx_g.number_of_nodes() < n:
        nx_g.add_node(_Point(random.random(), random.random(), random.random()))
    nodes = list(nx_g.nodes())
    for i in range(random.randrange(n, 20 * n)):
        u, v = random.sample(nodes, 2)
        nx_g.add_edge(u, v, weight=_distance(u, v))

    g = DistanceOracle(nx_g)
    for i in range(1000):
        s = random.choice(nodes)
        t = random.choice(nodes)
        path = g.find_shortest_path(s, t)
        nx_path = (
            networkx.dijkstra_path(nx_g, s, t)
            if networkx.has_path(nx_g, s, t)
            else None
        )
        assert path == nx_path, (_total_distance(path), _total_distance(nx_path))


if __name__ == "__main__":
    _test()
