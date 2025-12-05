def Newton_diagram(V):
    """
    return a list of all the faces of the Newton diagram
    """
    n = len(V[0])
    R = []
    for j in range(n):
        e = [0] * n
        e[j] = 1
        R.append(e)
    P = Polyhedron(vertices=V, rays=R)
    r = []
    for d in range(n):
        for f in P.faces(d):
            if len(f.rays()) == 0:
                r.append(f)
    return r

def sigma(f):
    n = len(f.vertices()[0])
    r = [False] * n
    for v in f.vertices():
        for j in range(n):
            if v[j] != 0:
                r[j] = True
    return r

def restricted_face(f):
    n = len(f.vertices()[0])
    s = sigma(f)
    V = [
            [ v[i] for i in range(n) if s[i] ]
        for v in f.vertices() ]
    return Polyhedron(V)

def Newton_number(V):
    G = Newton_diagram(V)
    n = len(V[0])
    r = 0;
    for f in G:
        s = sigma(f)
        d = len([ i for i in range(n) if s[i] ])
        r += (-1)**(n-d) * factorial(d) * Polyhedron(list(restricted_face(f).vertices()) + [[0]*d]).volume()
    return r + (-1)**n

def normv_red(f):
    fred = restricted_face(f)
    v = fred.vertices()[0]
    S,U,W = matrix([ vector(w) - vector(v) for w in fred.vertices()]).smith_form()
    return pos_int_prim(W.column(W.ncols()-1))

def normv(f):
    n = f.ambient_dimension()
    si = sigma(f)
    c = [ i for i in range(n) if si[i] ]
    d = [0] * n
    for i in range(len(c)):
        d[c[i]] = i
    nvred = normv_red(f)
    r = [ oo ] * n
    for i in c:
        r[i] = nvred[d[i]]
    return r

def mF(f):
    fred = restricted_face(f)
    v = fred.vertices()[0]
    return normv_red(fred).inner_product(vector(v))

def nF(f):
    fred = restricted_face(f)
    v = fred.vertices()[0]
    return min(normv_red(fred))

def maxax(f):
    return mF(f) / nF(f)

def is_coordinate(f):
    fred = restricted_face(f)
    return f.dimension() == fred.ambient_dimension() - 1

def alt_jacobian_polygon(V, d=-1):
    n = len(V[0]) - 1
    if d == -1:
        d = n
        G = Newton_diagram(V)
        r = [];
        for f in G:
            si = sigma(f)
            s = len([ i for i in range(n+1) if si[i] ]) - 1
            fred = restricted_face(f)
            v = fred.vertices()[0]
            S,U,W = matrix([ vector(w) - vector(v) for w in fred.vertices()]).smith_form()
            if S.ncols() - S.rank() == 1:
                nv = pos_int_prim(W.column(W.ncols()-1))
                mf = nv.inner_product(vector(v))
                ng = min(nv)
                vo = factorial(s+1) * Polyhedron(list(fred.vertices()) + [[0]*(s+1)]).volume() / mf
                r.append( [ (-1)**(n-s) * vo, mf, ng ] )
        return normalize_polyhedron(r)
    r = []
    for v in NCF(times_gen_lin(V)):
        s = len( [ i for i in range(n+1) if v[i] != oo ]) - 1
        mv = wt(v,V)
        nv = min(v)
        w = W_function(d,v,V)
        r.append( [(-1)**(n-s) * w, mv, nv] )
    return normalize_polyhedron(r)

def jacobian_polygon(V, d=-1):
    n = len(V[0]) - 1
    if d == -1:
        d = n
    return normalize_polyhedron( alt_jacobian_polygon(V, d) + alt_jacobian_polygon(V, d-1))

def pos_int_prim(v):
    if len(v) and v[0] < 0:
        return -v / gcd(v)
    else:
        return v/gcd(v)

def mod_slope(f):
    return f[1] / (f[1] + f[2])

def normalize_polyhedron(p):
    ms = list({ mod_slope(f) for f in p })
    ms.sort()
    r =  [ [ sum([ f[0] * gcd(f[1], f[2]) for f in p if mod_slope(f) == s ]),
             s.numerator(),
             s.denominator() - s.numerator() ]
           for s in ms if s != 1]
    return [ f for f in r if f[0] != 0 ]

def vpoly_length(p):
    return sum([ p[i][0] * p[i][1] for i in range(len(p)) ])

def vpoly_height(p):
    return sum([ p[i][0] * p[i][2] for i in range(len(p)) ])

def suspension_exp(V,N):
    n = len(V[0])
    return [ v + [0] for v in V ] + [ [0]*n + [N] ]

def suspension_type_two(V, iterate=1):
    if iterate == 0:
        return V
    if len(V) == 0:
        n = 0
    else:
        n = len(V[0])
    return suspension_type_two([ v + [0,0] for v in V ] + [ [0]*n + [1,1] ], iterate-1)

def Cap(S):
    l = len(S)
    if l == 0:
        return 1
    n = len(S[0])
    P = Polyhedron([vector([0]*n)] + [ sum( [vector(S[i]) for i in I] ) for I in subsets(range(l)) if len(I) > 0 ])
    return len([ p for p in P.integral_points() if P.relative_interior_contains(p) ])

def vol_by_Cap(S):
    l = len(S)
    if l == 0:
        return 1
    n = len(S[0])
    return sum([ Cap([ S[i] for i in I]) for I in subsets(range(l)) ])

def np_deg(N):
    NN = normalize_polyhedron(N)
    if len(NN) == 0:
        return -1
    f = NN[-1]
    if f[2] == 0:
        return "infinity"
    return f[1] / f[2]

def PC(V):
    G = Newton_diagram(V)
    return PolyhedralComplex([Polyhedron(f) for f in G])

def comb_newt_polyhedron(A, rel=False):
    n = A.ambient_dimension()
    r = []
    po = A.face_poset()
    for f in A.cell_iterator():
        if is_coordinate(f) and (rel==False or po.compare_elements(f,rel) in {0,1}):
            r.append([ (-1)**(n-f.dim()-1)/mF(f) , mF(f), nF(f) ])
    return normalize_polyhedron(r)

def max_deg_cnp(A):
    return max([ np_deg(comb_newt_polyhedron(A,f)) for f in A.cell_iterator() if Cap(f.vertices()) ] + [np_deg(comb_newt_polyhedron(A))])

def the_conjecture_tester(V, print_nonzeros=False):
    if not isolated_by_Kouchnirenko(V):
        print("Not isolated!")
        return True
    
    A = PC(V)
    B = A.subdivide(make_simplicial=True)

    maxdeg_facets = -1
    maxdeg_simplex = max_deg_cnp(B)
    loj_exp = Lojasiewicz_exponent(V)
    ajp = alt_jacobian_polygon(V)
    
    for f in A.cell_iterator():
        if is_coordinate(f):
            for g in B.cell_iterator():
                if is_contained_in(g,f) and Cap(g.vertices()) > 0 and np_deg(comb_newt_polyhedron(B,g)) >= maxax(f):
                    print("Example:")
                    print(f.vertices())
                    print(g.vertices())
                    maxdeg_facets = max( maxdeg_facets, maxax(f) )
  
    print("Lojasiewicz exponent:     ", loj_exp)
    print("maxdeg_facets minus one:  ", maxdeg_facets - 1)
    print("maxdeg_simplex minus one: ", maxdeg_simplex - 1)

    even_morse_point = False
    if ajp == []:
        even_morse_point = True
        print("Morse point in even number of variables")

    if np_deg(ajp) == maxdeg_facets or (even_morse_point and maxdeg_facets <= 2):
        print("The conjecture holds for this one.")
        return True
    else:
        print("Eureka, it is a counterexample!")
        return False

def isolated_by_Kouchnirenko(V):
    n = len(V[0])
    M = []
    for j in range(n):
        M.append([])
        for i in range(len(V)):
            e = [V[i][k] for k in range(n)]
            if e[j] != 0:
                e[j] = e[j] - 1
                M[j].append(e)
    for I in subsets(range(n)):
        c = 0
        for j in range(n):
            if len([x for x in M[j] if sum([ x[k] for k in I ]) == 0 ]) > 0:
                c += 1
        if c < n-len(I):
            return False
    return True

def times_gen_lin(V):
    if len(V) == 0:
        return []
    n = len(V[0])
    E = [ vector([0]*n) for i in range(n) ]
    for i in range(n):
        E[i][i] = 1
    return [ list(vector(v) + e) for v in V for e in E ]

def NCF(V):
    G = Newton_diagram(V)
    return [ normv(f) for f in G if is_coordinate(f) ]

def vector_pairing(u,v):
    n = len(u)
    r = 0
    for i in range(n):
        if u[i] == oo and v[i] > 0:
            return oo
        if v[i] == oo and u[i] > 0:
            return oo
        if u[i] != oo and v[i] != oo:
            r += u[i] * v[i]
    return r

def wt(v,V):
    return min([ vector_pairing(v,V[i]) for i in range(len(V)) ])

def W_function(d,v,V):
    n = len(v)-1
    A = K(v,V)
    s = len([ i for i in range(n+1) if v[i] != oo ]) - 1
    E = [ vector([0]*(n+1)) for i in range(n+1) ]
    for i in range(n+1):
        E[i][i] = 1
    B = K(v,E)
    c = n-d
    if d == n:
        return mv_factor_gen( [ A for i in range(s) ] )
    else:
        return sum([ binomial(k-1,c-1) * mv_factor_gen( [ A for i in range(s-k) ] + [ B for i in range(k) ]) for k in range(c,s+1) ])

def K(v, V):
    w = wt(v,V)
    vertices = [ a for a in V if vector_pairing(v,a) == w ]
    return Polyhedron(vertices)

def bin(k,n):
    b = []
    for i in range(0,n):
        b.append(k % 2)
        k = k // 2
    return b

def mv_factor(K):
    n = len(K)
    if n == 0:
        return 1
    v = 0
    for k in range(0,2**n):
        b = bin(k,n)
        N = Polyhedron([ [0]*n ])
        for j in range(0,n):
            if b[j]:
                N = N.minkowski_sum(K[j])
        v += (-1)**(n - sum(b)) * N.volume()
    return v

def mv_factor_gen(K):
    s = len(K)
    if s == 0:
        return 1
    n = K[0].ambient_dimension()
    M = FreeModule(ZZ,n)
    vertices = [ k.vertices() for k in K ]
    vertex_vectors = [ [ w.vector() for w in v ] for v in vertices ]
    differences = [ [ a - b[0] for a in b if a != b[0] ] for b in vertex_vectors ]
    all_diff = [ list(d) for v in differences for d in v ]
    N = M.submodule( all_diff ).saturation()
    if N.rank() != s:
        return 0
    coords = [ [vector([0]*s)] +  [ N.coordinate_vector(d) for d in v ] for v in differences ]
    polys = [ Polyhedron(c) for c in coords ]
    return mv_factor(polys)

def Lojasiewicz_exponent(V):
    if isolated_by_Kouchnirenko(V):
        return max(1, np_deg(alt_jacobian_polygon(V)) - 1)
    else:
        return oo

def ell_vector(B):
    d = B.dimension() + 1
    r = vector( [0]*(d+1) )
    r[d] = (-1)**(d-1)
    for f in B.cell_iterator():
        s = len( [ i for i in range(n+1) if sigma(f)[i] ] ) - 1
        e = excess(f)
        noG = f.dim() + 1
        for j in range(e+1):
            r[d-j] += (-1)**(d - noG + j - 1) * binomial(e,j)
    return r

def excess(f):
    n = f.ambient_dimension() - 1
    s = len( [ i for i in range(n+1) if sigma(f)[i] ] ) - 1
    return s - f.dimension()

def is_contained_in(G,F):
    r = True
    for v in G.vertices():
        if not F.contains(v):
            r = False
    return r
