# Successive Vertex Ordering Enumerator

This repository contains a Python implementation of the inclusion-exclusion method used to count successive vertex orderings of a finite connected simple graph.

The code follows the construction from the paper **"Successive vertex orderings of graphs"**.

## What the program computes

For a finite connected simple graph `G` on `n` vertices, the number of successive vertex orderings is computed from

`sigma(G) / n! = sum over independent sets I of (-1)^{|I|} * a(I)/n * b(I)`

where:

- `a(I) = n - |N[I]|`
- `b(∅) = 1`
- `b(I) = (1 / |N[I]|) * sum_{v in I} b(I \ {v})` for nonempty independent sets

The implementation builds the independent sets in lexicographic order, computes `b(I)` recursively, then sums the signed contributions.

## Files

- `svo_enumerator.py` — the main implementation

## Usage

```python
from svo_enumerator import successive_vertex_orderings

A = [
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0],
    [0, 1, 0, 1, 0],
    [0, 1, 1, 0, 1],
    [1, 0, 0, 1, 0],
]

print(successive_vertex_orderings(A))
```

For the worked example in the paper, this returns `60`.

## Assumptions

The implementation assumes:

- the input is a **simple graph**;
- the adjacency matrix is square, symmetric, and has zero diagonal;
- the graph is connected.

These assumptions match the main theorem in the paper.

## Relation to the paper

The code is a direct implementation of the proof strategy in the paper:

1. enumerate all independent sets;
2. compute `a(I) = n - |N[I]|`;
3. evaluate `b(I)` by the recursion;
4. sum `(-1)^{|I|} a(I)/n * b(I)` over all independent sets;
5. multiply by `n!`.

The classes in the code correspond to the notation in the paper as follows:

- `Subset` represents a vertex subset `I`;
- `Graph.neighborhood(I)` computes `N(I)`;
- `Node.a` stores `a(I)`;
- `Node.b` stores `b(I)`;
- `Node.value` stores the signed contribution to `sigma(G)/n!`.

## Sanity check against the paper

The included self-test uses the 5-vertex graph from Appendix A of the paper. The program returns `60`, which matches the worked example.

## Limitations

The algorithm is exact but exponential in the worst case, because the number of independent sets can be exponential in the number of vertices. That is exactly the complexity profile discussed in the paper.
