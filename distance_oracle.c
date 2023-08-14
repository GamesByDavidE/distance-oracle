#include "distance_oracle.h"
#include <math.h>
#include <stdlib.h>

enum { BUCKET_COUNT = 10000 };

struct Entry {
  struct Entry *previous;
  struct Entry *next;
};

typedef uint32_t Epoch;

struct Location {
  Coordinate x;
  Coordinate y;
  Coordinate z;
};

struct Node {
  struct Entry entry;
  struct Location location;
  Weight distance;
  NodeIndex parent;
  Epoch epoch;
  ArcIndex first_outgoing_arc;
};

struct Arc {
  Weight weight;
  NodeIndex head;
};

typedef uint32_t BucketIndex;

struct DistanceOracle {
  struct Node *nodes;
  struct Arc *arcs;
  struct Entry *buckets;
  Epoch epoch;
};

static void InitializeEntry(struct Entry *entry) {
  entry->previous = entry;
  entry->next = entry;
}

static void Splice(struct Entry *a, struct Entry *b) {
  struct Entry *temp = a->next;
  a->next = b->next;
  b->next = temp;
  a->next->previous = a;
  b->next->previous = b;
}

static void Relax(struct Node *node, Weight distance, NodeIndex parent,
                  Epoch epoch, BucketIndex *bucket_count,
                  struct Entry *buckets) {
  if (!(node->epoch < epoch || distance < node->distance))
    return;
  Splice(&node->entry, node->entry.previous);
  node->distance = distance;
  node->parent = parent;
  node->epoch = epoch;
  BucketIndex i = (BucketIndex)fmax(0, fmin(node->distance, BUCKET_COUNT - 1));
  if (*bucket_count <= i)
    *bucket_count = i + 1;
  Splice(&node->entry, &buckets[i]);
}

static bool IsEmpty(struct Entry *entry) { return entry == entry->previous; }

static struct Entry *Pop(struct Entry *entry) {
  struct Entry *front = entry->next;
  Splice(front, entry);
  return front;
}

static Weight Distance(struct Location a, struct Location b) {
  return fmax(fabs(a.x - b.x), fmax(fabs(a.y - b.y), fabs(a.z - b.z)));
}

#define NEW(pointer, count)                                                    \
  do {                                                                         \
    pointer = calloc(count, sizeof *pointer);                                  \
    if (pointer == NULL && count == 0)                                         \
      abort();                                                                 \
  } while (false)

void NewDistanceOracle(DistanceOracleHandle *graph_handle, NodeIndex node_count,
                       const ArcIndex first_outgoing_arcs[],
                       const NodeIndex arc_heads[]) {
  struct DistanceOracle *graph;
  NEW(graph, 1);

  NEW(graph->nodes, node_count + 1);
  for (NodeIndex node = 0; node < node_count; node++) {
    struct Node *n = &graph->nodes[node];
    InitializeEntry(&n->entry);
    n->first_outgoing_arc = first_outgoing_arcs[node];
  }
  graph->nodes[node_count].first_outgoing_arc = first_outgoing_arcs[node_count];

  ArcIndex arc_count = first_outgoing_arcs[node_count];
  NEW(graph->arcs, arc_count);
  for (ArcIndex arc = 0; arc < arc_count; arc++) {
    struct Arc *a = &graph->arcs[arc];
    a->weight = HUGE_VAL;
    a->head = arc_heads[arc];
  }

  NEW(graph->buckets, BUCKET_COUNT);
  for (BucketIndex i = 0; i < BUCKET_COUNT; i++)
    InitializeEntry(&graph->buckets[i]);

  *graph_handle = graph;
}

#undef NEW

void PlaceNode(DistanceOracleHandle graph_handle, NodeIndex node, Coordinate x,
               Coordinate y, Coordinate z) {
  struct DistanceOracle *graph = graph_handle;
  struct Location *location = &graph->nodes[node].location;
  location->x = x;
  location->y = y;
  location->z = z;
}

void LightenArc(DistanceOracleHandle graph_handle, ArcIndex arc,
                Weight new_weight) {
  struct DistanceOracle *graph = graph_handle;
  graph->arcs[arc].weight = new_weight;
}

bool FindShortestPath(DistanceOracleHandle graph_handle, NodeIndex s,
                      NodeIndex t, NodeIndex *path_node_count,
                      NodeIndex path_nodes[]) {
  struct DistanceOracle *graph = graph_handle;
  graph->epoch++;
  Epoch epoch = graph->epoch;
  BucketIndex bucket_count = 0;
  Relax(&graph->nodes[s], 0, s, epoch, &bucket_count, graph->buckets);
  BucketIndex front_index = 0;
  bool found_t = false;
  struct Location t_location = graph->nodes[t].location;
  while (true) {
    struct Entry *front = &graph->buckets[front_index];
    if (IsEmpty(front)) {
      if (found_t)
        break;
      front_index++;
      if (bucket_count <= front_index)
        return false;
      continue;
    }
    struct Node *n = (struct Node *)Pop(front);
    if (n->epoch < epoch)
      continue;
    NodeIndex node = n - graph->nodes;
    if (node == t)
      found_t = true;
    Weight distance_n_t = Distance(n->location, t_location);
    for (ArcIndex arc = n->first_outgoing_arc;
         arc < graph->nodes[node + 1].first_outgoing_arc; arc++) {
      struct Arc *a = &graph->arcs[arc];
      struct Node *h = &graph->nodes[a->head];
      Weight reduced_weight =
          a->weight - distance_n_t + Distance(h->location, t_location);
      Relax(h, n->distance + reduced_weight, node, epoch, &bucket_count,
            graph->buckets);
    }
  }

  NodeIndex k = 0;
  for (NodeIndex u = t; true; u = graph->nodes[u].parent) {
    path_nodes[k] = u;
    k++;
    if (u == s)
      break;
  }
  *path_node_count = k;
  for (NodeIndex i = 0, j = k - 1; i < j; i++, j--) {
    NodeIndex temp = path_nodes[i];
    path_nodes[i] = path_nodes[j];
    path_nodes[j] = temp;
  }
  return true;
}

#define DELETE(pointer)                                                        \
  do {                                                                         \
    free(pointer);                                                             \
    pointer = NULL;                                                            \
  } while (false)

void DeleteDistanceOracle(DistanceOracleHandle *graph_handle) {
  struct DistanceOracle *graph = *graph_handle;
  *graph_handle = NULL;
  DELETE(graph->buckets);
  DELETE(graph->arcs);
  DELETE(graph->nodes);
  DELETE(graph);
}

#undef DELETE
