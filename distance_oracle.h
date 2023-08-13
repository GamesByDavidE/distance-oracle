#include <stdbool.h>
#include <stdint.h>

typedef uint32_t NodeIndex;
typedef uint32_t ArcIndex;
typedef double Coordinate;
typedef double Weight;
typedef void *DistanceOracleHandle;

void NewDistanceOracle(DistanceOracleHandle *graph_handle, NodeIndex node_count,
                       const ArcIndex first_outgoing_arcs[],
                       const NodeIndex arc_heads[]);

void PlaceNode(DistanceOracleHandle graph_handle, NodeIndex node, Coordinate x,
               Coordinate y, Coordinate z);

void LightenArc(DistanceOracleHandle graph_handle, ArcIndex arc,
                Weight new_weight);

bool FindShortestPath(DistanceOracleHandle graph_handle, NodeIndex s,
                      NodeIndex t, NodeIndex *path_node_count,
                      NodeIndex path_nodes[]);

void DeleteDistanceOracle(DistanceOracleHandle *graph_handle);
