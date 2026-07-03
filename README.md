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

### `Polyomino`

Represents a graph using its adjacency matrix.

Provides methods for

- neighbourhood computation,
- independent set testing,
- graph operations used by the inclusion-exclusion algorithm.

Although the class is named `Polyomino`, it accepts any graph represented by an adjacency matrix.

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
graph = rectangular_grid_adj_matrix(3,2).tolist()

lattice = Lattice(Polyomino(graph))

coefficients = lattice.get_polynomial_coefficients()

print(coefficients)
```

---

## Computing the \(A_k\) Coefficients

```python
coeffs = lattice.get_polynomial_coefficients()

A_k = []

for k in range(len(coeffs)):
    value = kth_derivative_at_minus_one(coeffs, k)

    A_k.append(value)

print(A_k)
```

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

The input graph should satisfy

- the adjacency matrix is square;
- entries are either 0 or 1;
- the graph is simple (no self-loops);
- the graph is undirected (symmetric adjacency matrix).

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
