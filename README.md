# mesh_lod232

Mesh lightweighting (surface simplification) for digital-workpiece assets,
implemented from scratch in C++20 on top of Eigen 3.4 (header-only, expected
at `third_party/eigen3`). Namespace: `mesh_lod232`.

## Algorithm

Quadric-error metric edge collapse (Garland–Heckbert style):

- Each vertex starts with the sum of the plane quadrics of the **initial**
  faces (unit-normal planes). Quadrics are merged additively on collapse and
  are never re-initialised from the current surface.
- Collapse candidates are restricted to current mesh edges. The optimal
  position minimises the summed quadric; if the linear system is singular,
  the two endpoints and the midpoint are compared instead.
- Candidates are processed by increasing error; ties are broken by the
  stable original vertex ids. Candidate costs and legality always reflect
  the adjacency after the previous collapse (lazy invalidation via
  per-vertex version counters; stale heap entries are discarded).

## Topology and safety rules

- All **initial boundary vertices are locked**: never moved, never removed.
- An edge collapse is accepted only if it satisfies the topological **link
  condition**, removes exactly the faces adjacent to the edge, and updates
  the remaining neighbouring faces.
- Collapses that would create duplicate faces, zero-area faces, or normal
  flips are rejected, as are collapses that would delete an entire connected
  component. Different components are never connected (collapses happen
  along existing edges only).
- No global self-intersection guard and no optimality guarantee, by design.

## Input validation

`simplify` throws `mesh_lod232::ValidationError` on: non-finite coordinates,
out-of-range indices, duplicate faces, zero-area faces, inconsistent
winding, non-manifold edges, non-manifold vertices, or a zero target face
count. Multiple connected components and open boundaries are allowed. The
input containers are never modified. If the target is not smaller than the
current face count, the input is returned unchanged.

## Termination

Iteration stops when the face count is at most the target
(`kTargetReached`), or earlier when no legal candidate remains
(`kNoLegalCandidate`; the current mesh is kept and the result states that
the target was not reached). The result reports compacted vertices,
remapped triangle indices, the actual face count, the collapse count, and
the stop reason. Winding and locked boundary coordinates are preserved.

## Layout

- `include/mesh_lod232/mesh_lod.h` — public API (`simplify`, `SimplifyResult`).
- `src/validation.cpp` — input validation.
- `src/mesh.cpp` — mesh topology: adjacency, link condition, edge collapse,
  boundary (locked) vertices.
- `src/quadric.cpp` — quadric error metrics and optimal placement.
- `src/simplify.cpp` — iterative collapse driver.
- `examples/open_surface.cpp` — open-surface demo (curved grid with a hole
  plus a second component).
- `tests/selftest.cpp` — self tests.

## Build, test, example

```sh
make            # builds build/libmesh_lod232.a, example and tests
make test       # runs the self tests
make run-example
```

Example output:

```
input : 451 vertices, 776 faces (2 connected components, open boundaries)
output: 114 vertices, 120 faces
collapses: 328
stop reason: target face count reached
boundary vertices preserved: 106 / 106 (coordinates unchanged)
```

## Usage

```cpp
#include "mesh_lod232/mesh_lod.h"

std::vector<mesh_lod232::Vec3> vertices = /* ... */;
std::vector<mesh_lod232::Tri>  faces    = /* ... */;
auto result = mesh_lod232::simplify(vertices, faces, /*target_faces=*/1000);
```
