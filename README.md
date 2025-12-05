# jacobian_polygon

This script implements some of the formulas in arXiv:2509.06150 for
Newton nondegenerate singlarities. The function
`alt_jacobian_polyhedron` returns the alternating Jacobian polygon in
the form of a list of triples. Each triple `[a,m,n]` represents a term
`a {m//n}`.

As an example, let us consider the singularity `f(x,y,z) = x^2 + y^3 +
z^5`. As input, we give a list containing the vertices of the Newton
diagram.

```
sage: load("jacobian_polygon.sage")
sage: V = [ [2,0,0], [0,3,0], [0,0,5] ]
sage: alt_jacobian_polygon(V)
[[1, 2, 1], [-1, 3, 1], [2, 5, 1]]
```

As a result, the alternating Jacobian polygon is `{2//1} - {3//1} + 2
{5//1}`.
