"""
Title: Enumerate successive vertex orderings of a finite simple graph.

1. This script implements the inclusion--exclusion method described in the paper "Successive vertex orderings of graphs".
2. Calculate successive ordering polynomial of the graph.
3. Calculate Ak values by taking derivative of the polynomial.
4. Illustrate these calculations on example graphs.
5. Apply to jigsaw puzzles.

Core identity used by the implementation:
    sigma(G) / n! = sum_{I independent} (-1)^{|I|} * a(I)/n * b(I)
where
    a(I) = n - |N[I]|,
    b(∅) = 1,
    b(I) = (1 / |N[I]|) * sum_{v in I} b(I minus {v}) for I != ∅.

Author: Hadi, Prarthana
Date: 09 Feb 2026
Version: 1.0
"""

#========================================== Main script ==========================================#
# Import packages
from __future__ import annotations
from fractions import Fraction
import math
import numpy as np

class Subset:
    vector: list[int]
    n: int

    def __init__(self, vector: list[int]):
        self.vector = vector
        for i in vector:
            if i not in (0, 1):
                raise ValueError("indicator vector entries must be 0 or 1")
        self.n = len(vector)

    def add_element(self, pos: int):
        """Add element to pos."""
        if not 0 <= pos < self.n:
            raise IndexError(f"vertex position {pos} out of range for subset of size {self.n}")
        self.vector[pos] = 1

    def remove_element(self, pos: int):
        """Remove element at pos."""
        if not 0 <= pos < self.n:
            raise IndexError(f"vertex position {pos} out of range for subset of size {self.n}")
        self.vector[pos] = 0

    def __add__(self, other: Subset) -> Subset:
        """Gives the union of two subsets."""
        if self.n != other.n:
            raise ValueError(
                f"cannot take union of subsets over ground sets of different sizes ({self.n} vs {other.n})"
            )
        result_list = [self.vector[i] | other.vector[i] for i in range(self.n)]
        return Subset(result_list)

    def __mul__(self, other: Subset) -> Subset:
        """Gives the intersection of two subsets."""
        if self.n != other.n:
            raise ValueError(
                f"cannot intersect subsets over ground sets of different sizes ({self.n} vs {other.n})"
            )
        result_list = [self.vector[i] & other.vector[i] for i in range(self.n)]
        return Subset(result_list)

    def __abs__(self) -> int:
        """Returns size of subset."""
        return sum(self.vector)

    def get_children(self) -> list[Subset]:
        """Get all subsets with one more element."""
        children = []
        for i in range(self.n):
            if self.vector[i] == 0:
                child = Subset(self.vector.copy())
                child.add_element(i)
                children.append(child)
        return children

    def get_parents(self) -> list[Subset]:
        """Get all subsets with one less element."""
        parents = []
        for i in range(self.n):
            if self.vector[i] == 1:
                parent = Subset(self.vector.copy())
                parent.remove_element(i)
                parents.append(parent)
        return parents

    def __getitem__(self, i: int):
        return self.vector[i]

    def __setitem__(self, i: int, x: int):
        if x not in (0, 1):
            raise ValueError("indicator vector entries must be 0 or 1")
        self.vector[i] = x

    @property
    def largest_element(self) -> int:
        if abs(self) == 0:
            return -1
        else:
            return max([i for i in range(self.n) if self.vector[i] == 1])

    def get_next_subsets(self) -> list[Subset]:
        """Get all children subsets with one element added to right."""
        children = []
        largest_element = self.largest_element
        for i in range(largest_element + 1, self.n):
            child = Subset(self.vector.copy())
            child.add_element(i)
            children.append(child)
        return children

    def __eq__(self, other) -> bool:
        return (tuple(self.vector) == tuple(other.vector))

    def __hash__(self) -> int:
        return hash(tuple(self.vector))

    def __str__(self) -> str:
        string = ''
        for i in range(self.n):
            string += str(self.vector[i])
        return string

    def __repr__(self) -> str:
        return str(self)

def make_empty_subset(n: int) -> Subset:
    return Subset([0 for i in range(n)])

class SimpleGraph:
    adj_matrix: list[list[int]]
    n: int

    def __init__(self, adj_matrix: list[list[int]]):
        self.n = len(adj_matrix)
        for i in range(self.n):
            if len(adj_matrix[i]) != self.n:
                raise ValueError("adjacency matrix must be square")
            if adj_matrix[i][i] != 0:
                raise ValueError(f"self-loop at vertex {i}: the graph must be simple")
            for j in range(self.n):
                if adj_matrix[i][j] not in (0, 1):
                    raise ValueError("adjacency matrix entries must be 0 or 1")
                if adj_matrix[i][j] != adj_matrix[j][i]:
                    raise ValueError(
                        f"adjacency matrix must be symmetric: entries ({i},{j}) and ({j},{i}) differ"
                    )
        self.adj_matrix = adj_matrix

    def neighborhood(self, subset: Subset) -> Subset:
        neighborhood = make_empty_subset(self.n)
        for i in range(self.n):
            if subset[i] == 1:
                neighborhood += Subset(self.adj_matrix[i])
        return neighborhood

    def __mul__(self, subset: Subset) -> Subset:
        return self.neighborhood(subset)

    def is_independent(self, subset: Subset) -> bool:
        return abs((self * subset) * subset) == 0

    def __xor__(self, subset: Subset) -> bool:
        return self.is_independent(subset)

class Node:
    independent_set: Subset
    neighborhood: Subset
    n: int # number of vertices of the graph
    a: int # |V| - |N(I) + I|
    b: Fraction | None # parameter to be computed
    d: int # subset size
    children: list[Node]
    parents: list[Node]
    is_root: bool

    def __init__(self, subset: Subset, neighborhood: Subset):
        # Internal invariants (deliberately asserts, not exceptions): Node is only
        # constructed by Lattice, which guarantees both conditions.
        assert subset.n == neighborhood.n  # ground sets agree by construction
        assert abs(subset * neighborhood) == 0  # I and N(I) disjoint, i.e. I is independent
        self.independent_set = subset
        self.neighborhood = neighborhood
        self.n = self.independent_set.n
        self.a = self.n - abs(neighborhood) - abs(self.independent_set)
        self.d = abs(self.independent_set)
        self.children = []
        self.parents = []
        if abs(self.independent_set) == 0:
            self.is_root = True
            self.b = Fraction(1)
        else:
            self.is_root = False
            self.b = None # to be calculated

    def __eq__(self, other) -> bool:
        return self.independent_set == other.independent_set

    def __hash__(self) -> int:
        return hash(self.independent_set)

    @property
    def value(self) -> Fraction:
        # Internal invariant: Lattice computes b level-by-level before value is read.
        assert self.b is not None, "internal invariant: b must be computed before value is read"
        return Fraction(self.a) * Fraction(self.b) * Fraction(1, self.n) * (-1) ** self.d

    def calculate_b(self) -> Fraction:
        self.b = Fraction(0)
        for parent in self.parents:
            self.b += parent.b / Fraction(self.n - self.a)
        return self.b
    
class Lattice:
    n: int # number of vertices of the graph
    graph: SimpleGraph
    root: Node
    levels: dict[int, list[Node]] # returns nodes with that many elements
    nodes: dict[int, Node] # key is hash value of independent set tuple

    def __init__(self, graph: SimpleGraph):
        self.n = graph.n
        self.graph = graph
        self.root = Node(make_empty_subset(self.n), make_empty_subset(self.n))
        self.levels = {0: [self.root]}
        self.nodes = {make_empty_subset(self.n): self.root}
        for d in range(1, self.n + 1):  # up to size n: for connected graphs alpha(G) <= n-1, but edgeless inputs have an independent set of size n
            self.levels[d] = []
            for node in self.levels[d - 1]:
                next_subsets = node.independent_set.get_next_subsets()
                for subset in next_subsets:
                    if not self.graph ^ subset:
                        continue
                    # N(I + {v}) = N(I) + N(v): reuse the parent's stored neighborhood instead of recomputing it
                    new_neighborhood = node.neighborhood + Subset(self.graph.adj_matrix[subset.largest_element])
                    new_node = Node(subset, new_neighborhood)
                    self.levels[d].append(new_node)
                    self.nodes[new_node.independent_set] = new_node
                    new_node.parents = [self.nodes[parent] for parent in new_node.independent_set.get_parents()]
                    for parent_node in new_node.parents:
                        parent_node.children.append(new_node)
                    new_node.calculate_b()
            if self.levels[d] == []:
                self.levels.pop(d)
                break

    def sum_all_values(self) -> Fraction:
        total_sum = Fraction(0)
        for node in self.nodes.values():
            total_sum += node.value
        return total_sum

    def get_polynomial_coefficients(self, a_included = True) -> list[Fraction]:
        coefficients = [Fraction(0) for i in range(self.n + 1)]
        for node in self.nodes.values():
            if a_included:
                coefficients[node.d] += node.value * (-1)**node.d
            else:
                coefficients[node.d] += node.b
        return coefficients   
        
#--------------------------------- SVO calculation -------------------------------------#
def successive_vertex_orderings(adj_matrix):
    """
    Given a graph (adjacency matrix), count the number of distinct growth sequences (successive vertex orderings).
    """

    graph = SimpleGraph(adj_matrix)
    lattice = Lattice(graph)
    fraction = lattice.sum_all_values()

    result = fraction * math.factorial(graph.n)
    if result.denominator != 1:
        raise ArithmeticError(
            f"internal error: sigma(G) evaluated to the non-integer {result}"
        )
    return result.numerator

#--------------------- Create adjacency matrix of an mxn puzzle ---------------------------#
def rectangular_grid_adj_matrix(m, n):
    """
    Returns adjacency matrix for an m x n grid (m rows, n columns) with 4-neighbour connectivity.
    """
    A = np.zeros((m * n, m * n), dtype=int)
    
    for i in range(m):
        for j in range(n):
            idx = i * n + j
            
            if i > 0: A[idx, (i - 1) * n + j] = 1 # Up
            if i < m - 1: A[idx, (i + 1) * n + j] = 1 # Down
            if j > 0: A[idx, i * n + (j - 1)] = 1 # Left
            if j < n - 1: A[idx, i * n + (j + 1)] = 1 # Right
                
    return A

#---------------------- Return derivative of the polynomial at x=-1 -------------------------# 
def kth_derivative_at_minus_one(coeffs, k):
    total = 0
    for n in range(k, len(coeffs)):
        falling = math.factorial(n) // math.factorial(n - k)
        total += coeffs[n] * falling * ((-1) ** (n - k))
    return total


if __name__ == "__main__":

    #===========================================================================================#
    #----------------------------------------- EXAMPLES ----------------------------------------#
    #===========================================================================================#

    #-------------------------------------------------------------------------#
    #  cycle graph with additional chord
    adj_matrix = [
        [0, 1, 0, 0, 1],  # Connections: 1-2, 1-5
        [1, 0, 1, 1, 0],  # Connections: 2-1, 2-3, 2-4
        [0, 1, 0, 1, 0],  # Connections: 3-2, 3-4
        [0, 1, 1, 0, 1],  # Connections: 4-2, 4-3, 4-5
        [1, 0, 0, 1, 0]   # Connections: 5-1, 5-4
    ]
    print(successive_vertex_orderings(adj_matrix))

    #----------------------------------------------------------------------------#
    # Example graph: 20 vertices
    adj_matrix = [
      [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 
      [1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
      [0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
      [1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
      [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
      [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
      [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
      [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
      [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], 
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], 
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0], 
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1], 
      [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    ]
    print(successive_vertex_orderings(adj_matrix))

    #-----------------------------------------------------------------------------#

    #======================================= Jigsaw puzzle =====================================#
    rows = 3; cols = 2
    size = rows*cols # 3x2 puzzle
    adj_matrix = rectangular_grid_adj_matrix(rows,cols).tolist() # convert to adjacency matrix

    # Build the independent-set lattice once; all queries below reuse it
    lattice = Lattice(SimpleGraph(adj_matrix))

    # Return successive vertex orderings
    svo = lattice.sum_all_values()
    print('SVO =', svo*math.factorial(size))

    # Return successive vertex ordering polynomial
    polynomial_coeffs = lattice.get_polynomial_coefficients(a_included = True)
    print(polynomial_coeffs)

    # Get Ak coefficients by taking derivative of svo polynomial
    A_k_list = []
    for k in range(size + 1):
        F_prime_k_at_minus_one = kth_derivative_at_minus_one(polynomial_coeffs, k)
        A_k = Fraction(math.factorial(size)) * F_prime_k_at_minus_one / math.factorial(k)
        if A_k.denominator != 1:
            raise ArithmeticError(
                f"internal error: A_{k} evaluated to the non-integer {A_k}"
            )
        A_k_list.append(A_k.numerator)
    print(A_k_list)

    # Consistency checks (Taylor expansion of the polynomial around x = -1):
    #   A_0 = sigma(G)   and   sum_k A_k = n!
    if A_k_list[0] != svo * math.factorial(size):
        raise ArithmeticError("consistency check failed: A_0 != sigma(G)")
    if sum(A_k_list) != math.factorial(size):
        raise ArithmeticError("consistency check failed: sum of A_k != n!")
    print("Consistency checks passed: A_0 = sigma(G) and sum(A_k) = n!")
