# Successive Vertex Ordering Enumerator

This repository contains a Python implementation of the inclusion-exclusion algorithm described in the paper **"Successive Vertex Orderings of Graphs"**.

The program computes the exact number of successive vertex orderings of a graph by enumerating all independent sets and evaluating the recursive quantities appearing in the inclusion-exclusion formula.

In addition, the implementation can

- compute the successive vertex ordering polynomial,
- compute the derivatives of the polynomial,
- calculate the corresponding \(A_k\) coefficients, and
- generate rectangular grid graphs for jigsaw puzzle applications.

---

## Mathematical Background

For a graph \(G=(V,E)\) with \(n=|V|\),

```
σ(G)/n! = Σ (-1)^|I| · (a(I)/n) · b(I)
```

where the summation is over all independent sets \(I\),

```
a(I) = n - |N[I]|
```

and

```
b(∅) = 1
```

```
b(I) = (1/|N[I]|) Σ b(I \ {v})
```

where the second summation is over all vertices \(v\in I\).

The implementation evaluates this identity exactly using rational arithmetic (`fractions.Fraction`).

---

## Features

The program provides:

- exact computation of successive vertex orderings;
- enumeration of all independent sets;
- construction of the independent-set lattice;
- computation of the recursive values \(b(I)\);
- computation of the successive vertex ordering polynomial;
- computation of polynomial derivatives evaluated at \(x=-1\);
- computation of the \(A_k\) coefficients;
- generation of rectangular grid graphs.

---

## Repository Structure

```
.
├── svo_enumerator.py
└── README.md
```

---

## Main Classes

### `Subset`

Represents a subset of vertices using a binary indicator vector.

Provides methods for

- adding and removing vertices,
- computing unions and intersections,
- generating parent and child subsets,
- lexicographic subset generation.

---

### `SimpleGraph`

Represents a finite simple undirected graph using its adjacency matrix.

Provides methods for

- neighbourhood computation,
- independent set testing,
- graph operations used by the inclusion-exclusion algorithm.

The constructor validates its input and raises a `ValueError` if the matrix is not square, not symmetric, contains entries other than 0 and 1, or has a nonzero diagonal (self-loop).

---

### `Node`

Represents one independent set in the lattice.

Each node stores

- the independent set,
- its neighbourhood,
- the quantity

```
a(I)=n-|N[I]|
```

- the recursively computed value

```
b(I)
```

- the corresponding contribution to the inclusion-exclusion sum.

---

### `Lattice`

Constructs the lattice of independent sets.

It computes

- all independent sets,
- parent-child relationships,
- recursive \(b(I)\) values,
- the inclusion-exclusion sum,
- successive vertex ordering polynomial coefficients.

---

## Functions

### `successive_vertex_orderings(adj_matrix)`

Returns the exact number of successive vertex orderings of the graph.

Example

```python
from svo_enumerator import successive_vertex_orderings

A = [
    [0,1,0,0,1],
    [1,0,1,1,0],
    [0,1,0,1,0],
    [0,1,1,0,1],
    [1,0,0,1,0]
]

print(successive_vertex_orderings(A))
```

Output

```
60
```

---

### `rectangular_grid_adj_matrix(m, n)`

Constructs the adjacency matrix of an \(m\times n\) rectangular grid using 4-neighbour connectivity.

Example

```python
from svo_enumerator import rectangular_grid_adj_matrix

graph = rectangular_grid_adj_matrix(3,2)
```

---

### `kth_derivative_at_minus_one(coeffs, k)`

Evaluates the \(k\)-th derivative of the successive vertex ordering polynomial at

```
x=-1.
```

This is used to compute the \(A_k\) coefficients.

---

## Computing the Successive Vertex Ordering Polynomial

```python
from svo_enumerator import SimpleGraph, Lattice, rectangular_grid_adj_matrix

adj_matrix = rectangular_grid_adj_matrix(3,2).tolist()

lattice = Lattice(SimpleGraph(adj_matrix))

coefficients = lattice.get_polynomial_coefficients()

print(coefficients)
```

---

## Computing the \(A_k\) Coefficients

The coefficient \(A_k\) counts the orderings in which exactly \(k\) non-first
vertices appear before all of their neighbours, and is given by

```
A_k = (n!/k!) · P^(k)(-1)
```

where \(P\) is the successive vertex ordering polynomial.

```python
import math
from svo_enumerator import kth_derivative_at_minus_one

# continuing from the previous example
coeffs = lattice.get_polynomial_coefficients()
n = len(coeffs) - 1

A_k = []
for k in range(len(coeffs)):
    value = math.factorial(n) * kth_derivative_at_minus_one(coeffs, k) / math.factorial(k)
    A_k.append(int(value))

print(A_k)
```

In particular \(A_0 = σ(G)\), and the values satisfy \(Σ_k A_k = n!\).

---

## Included Examples

The script includes examples demonstrating

- the five-vertex graph from the paper;
- a twenty-vertex graph;
- a \(3\times2\) rectangular grid graph;
- computation of successive vertex orderings;
- computation of the successive vertex ordering polynomial;
- computation of the \(A_k\) coefficients.

---

## Assumptions

The input graph must satisfy

- the adjacency matrix is square;
- entries are either 0 or 1;
- the graph is simple (no self-loops);
- the graph is undirected (symmetric adjacency matrix).

These conditions are checked by the `SimpleGraph` constructor, which raises a `ValueError` if any of them fails.

The implementation is intended for connected graphs, consistent with the main theorem in the paper.

---

## Complexity

The algorithm is exact.

Its running time is exponential in the worst case because it explicitly enumerates all independent sets of the graph. This matches the complexity discussed in the accompanying paper.

---

## Dependencies

The implementation requires

- Python 3.10+
- NumPy

Install NumPy using

```bash
pip install numpy
```

The Python standard library modules

- `math`
- `fractions`
- `typing`

are also used.

---

## Citation

If you use this implementation in academic work, please cite

> **Successive Vertex Orderings of Graphs**.
