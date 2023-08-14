extern "C" {
#include "distance_oracle.h"
}

#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

namespace {

struct Location {
  Coordinate x = 0;
  Coordinate y = 0;
  Coordinate z = 0;
};

Coordinate Distance(Location a, Location b) {
  return std::max({
      a.x - b.x,
      a.y - b.y,
      a.z - b.z,
  });
}

using Epoch = uint32_t;

struct ParentEntry {
  Weight distance = std::numeric_limits<Weight>::quiet_NaN();
  NodeIndex parent = std::numeric_limits<NodeIndex>::max();
  Epoch epoch = 0;
};

class ShortestPathTree {
public:
  Weight distance(NodeIndex child) const { return entries_[child].distance; }
  NodeIndex parent(NodeIndex child) const { return entries_[child].parent; }
  Epoch NewEpoch(NodeIndex node_count);
  bool UpdateParent(NodeIndex child, ParentEntry new_entry);

private:
  Epoch current_epoch_ = 0;
  std::vector<ParentEntry> entries_;
};

Epoch ShortestPathTree::NewEpoch(NodeIndex node_count) {
  if (current_epoch_ == std::numeric_limits<Epoch>::max()) {
    entries_.clear();
    current_epoch_ = 1;
  } else {
    ++current_epoch_;
  }
  entries_.resize(node_count);
  return current_epoch_;
}

bool ShortestPathTree::UpdateParent(NodeIndex child, ParentEntry new_entry) {
  ParentEntry &old_entry = entries_[child];
  if (!(old_entry.epoch < new_entry.epoch ||
        new_entry.distance < old_entry.distance))
    return false;
  old_entry = new_entry;
  return true;
}

class DistanceOracle {
public:
  explicit DistanceOracle(NodeIndex node_count,
                          std::vector<ArcIndex> first_outgoing_arcs,
                          std::vector<NodeIndex> arc_heads);
  void PlaceNode(NodeIndex node, Location location);
  void LightenArc(ArcIndex arc, Weight new_weight);
  std::optional<std::vector<NodeIndex>> FindShortestPath(NodeIndex s,
                                                         NodeIndex t);

private:
  using Entry = std::pair<Weight, NodeIndex>;
  Weight LowerBound(NodeIndex s, NodeIndex t) const;

  NodeIndex node_count_;
  std::vector<Location> locations_;
  std::vector<ArcIndex> first_outgoing_arcs_;
  std::vector<NodeIndex> arc_heads_;
  std::vector<Weight> arc_weights_;
  ShortestPathTree tree_;
  std::vector<std::vector<Entry>> buckets_;
};

DistanceOracle::DistanceOracle(NodeIndex node_count,
                               std::vector<ArcIndex> first_outgoing_arcs,
                               std::vector<NodeIndex> arc_heads)
    : node_count_{node_count}, locations_{node_count},
      first_outgoing_arcs_{std::move(first_outgoing_arcs)},
      arc_heads_{std::move(arc_heads)},
      arc_weights_(arc_heads_.size(), std::numeric_limits<Weight>::infinity()) {
  assert(!first_outgoing_arcs_.empty());
  assert(first_outgoing_arcs_.size() - 1 == node_count);
  assert(first_outgoing_arcs_.front() == 0);
  assert(
      std::is_sorted(first_outgoing_arcs_.begin(), first_outgoing_arcs_.end()));
  assert(first_outgoing_arcs_.back() == arc_heads_.size());
  for (NodeIndex node : arc_heads_)
    assert(node < node_count);
}

void DistanceOracle::PlaceNode(NodeIndex node, Location location) {
  assert(node < node_count_);
  locations_[node] = location;
}

void DistanceOracle::LightenArc(ArcIndex arc, Weight new_weight) {
  assert(arc < arc_heads_.size());
  assert(0 <= new_weight);
  Weight &old_weight = arc_weights_[arc];
  assert(new_weight <= old_weight);
  old_weight = new_weight;
}

std::optional<std::vector<NodeIndex>>
DistanceOracle::FindShortestPath(NodeIndex s, NodeIndex t) {
  assert(s < node_count_);
  assert(t < node_count_);
  Epoch epoch = tree_.NewEpoch(node_count_);
  {
    tree_.UpdateParent(s, {.distance = 0, .parent = s, .epoch = epoch});
    buckets_.resize(1);
    buckets_.front() = {{0, s}};
    long front = 0;
    bool found_t = false;
    for (;;) {
      if (buckets_[front].empty()) {
        if (found_t)
          break;
        ++front;
        if (buckets_.size() <= front)
          return std::nullopt;
        continue;
      }
      auto [dist_tail, tail] = buckets_[front].back();
      buckets_[front].pop_back();
      if (dist_tail != tree_.distance(tail))
        continue;
      if (tail == t)
        found_t = true;
      for (ArcIndex arc = first_outgoing_arcs_[tail];
           arc < first_outgoing_arcs_[tail + 1]; ++arc) {
        NodeIndex head = arc_heads_[arc];
        Weight dist_head =
            dist_tail +
            (arc_weights_[arc] - LowerBound(tail, t) + LowerBound(head, t));
        if (tree_.UpdateParent(
                head,
                {.distance = dist_head, .parent = tail, .epoch = epoch})) {
          auto i = static_cast<decltype(buckets_)::size_type>(dist_head);
          if (buckets_.size() <= i)
            buckets_.resize(i + 1);
          buckets_[i].push_back({dist_head, head});
        }
      }
    }
  }
  std::vector<NodeIndex> path = {t};
  while (path.back() != s)
    path.push_back(tree_.parent(path.back()));
  std::reverse(path.begin(), path.end());
  return std::move(path);
}

Weight DistanceOracle::LowerBound(NodeIndex s, NodeIndex t) const {
  assert(s < node_count_);
  assert(t < node_count_);
  return Distance(locations_[s], locations_[t]);
}

} // namespace

extern "C" {
void NewDistanceOracle(DistanceOracleHandle *graph_handle, NodeIndex node_count,
                       const ArcIndex first_outgoing_arcs[],
                       const NodeIndex arc_heads[]) {
  assert(graph_handle != nullptr);
  assert(*graph_handle == nullptr);
  assert(node_count < std::numeric_limits<NodeIndex>::max());
  assert(first_outgoing_arcs != nullptr);
  ArcIndex arc_count = first_outgoing_arcs[node_count];
  assert(arc_heads != nullptr || arc_count == 0);
  *graph_handle = new DistanceOracle{
      node_count,
      {first_outgoing_arcs, first_outgoing_arcs + node_count + 1},
      {arc_heads, arc_heads + arc_count},
  };
}

void PlaceNode(DistanceOracleHandle graph_handle, NodeIndex node, Coordinate x,
               Coordinate y, Coordinate z) {
  assert(graph_handle != nullptr);
  static_cast<DistanceOracle *>(graph_handle)->PlaceNode(node, {x, y, z});
}

void LightenArc(DistanceOracleHandle graph_handle, ArcIndex arc,
                Weight new_weight) {
  assert(graph_handle != nullptr);
  static_cast<DistanceOracle *>(graph_handle)->LightenArc(arc, new_weight);
}

bool FindShortestPath(DistanceOracleHandle graph_handle, NodeIndex s,
                      NodeIndex t, NodeIndex *path_node_count,
                      NodeIndex path_nodes[]) {
  assert(graph_handle != nullptr);
  assert(path_node_count != nullptr);
  assert(path_nodes != nullptr);
  auto path =
      static_cast<DistanceOracle *>(graph_handle)->FindShortestPath(s, t);
  if (!path)
    return false;
  std::copy(path->begin(), path->end(), path_nodes);
  *path_node_count = path->size();
  return true;
}

void DeleteDistanceOracle(DistanceOracleHandle *graph_handle) {
  assert(graph_handle != nullptr);
  delete static_cast<DistanceOracle *>(*graph_handle);
  *graph_handle = nullptr;
}
}
