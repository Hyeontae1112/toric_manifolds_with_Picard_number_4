#!/usr/bin/env python
# coding: utf-8

# In[83]:


from Simplicial_Complex import *
import sympy as sp
import z3


# In[213]:


from sympy.core import Mul, Expr, Add, Pow, Symbol, Number
def sympy_to_z3_real_var(sympy_var_list):
    z3_vars = []
    z3_var_map = {}
    
    for var in sympy_var_list:
        name = var.name
        z3_var = z3.Real(name)
        z3_var_map[name] = z3_var
        z3_vars.append(z3_var)

    return z3_vars, z3_var_map

def _sympy_to_z3_rec_real(var_map, e):
    'recursive call for sympy_to_z3()'

    rv = None

    if not isinstance(e, Expr):
        raise RuntimeError("Expected sympy Expr: " + repr(e))

    if isinstance(e, Symbol):
        rv = var_map.get(e.name)

        if rv == None:
            raise RuntimeError("No var was corresponds to symbol '" + str(e) + "'")

    elif isinstance(e, Number):
        rv = float(e)
    elif isinstance(e, Mul):
        rv = _sympy_to_z3_rec_real(var_map, e.args[0])

        for child in e.args[1:]:
            rv *= _sympy_to_z3_rec_real(var_map, child)
    elif isinstance(e, Add):
        rv = _sympy_to_z3_rec_real(var_map, e.args[0])

        for child in e.args[1:]:
            rv += _sympy_to_z3_rec_real(var_map, child)
    elif isinstance(e, Pow):
        term = _sympy_to_z3_rec_real(var_map, e.args[0])
        exponent = _sympy_to_z3_rec_real(var_map, e.args[1])

        if exponent == 0.5:
            # sqrt
            rv = z3.Sqrt(term)
        else:
            rv = term**exponent

    if rv == None:
        raise RuntimeError("Type '" + str(type(e)) + "' is not yet implemented for convertion to a z3 expresion. " + \
                            "Subexpression was '" + str(e) + "'.")

    return rv

def sympy_to_z3_eq_real(z3_var_map, sympy_exp):

    return _sympy_to_z3_rec_real(z3_var_map, sympy_exp)
    
def sympy_to_z3_var(sympy_var_list):
    z3_vars = []
    z3_var_map = {}
    
    for var in sympy_var_list:
        name = var.name
        z3_var = z3.Int(name)
        z3_var_map[name] = z3_var
        z3_vars.append(z3_var)

    return z3_vars, z3_var_map


def sympy_to_z3_eq(z3_var_map, sympy_exp):

    return _sympy_to_z3_rec(z3_var_map, sympy_exp)

def _sympy_to_z3_rec(var_map, e):
    'recursive call for sympy_to_z3()'

    rv = None

    if not isinstance(e, Expr):
        raise RuntimeError("Expected sympy Expr: " + repr(e))

    if isinstance(e, Symbol):
        rv = var_map.get(e.name)

        if rv == None:
            raise RuntimeError("No var was corresponds to symbol '" + str(e) + "'")

    elif isinstance(e, Number):
        rv = int(e)
    elif isinstance(e, Mul):
        rv = _sympy_to_z3_rec(var_map, e.args[0])

        for child in e.args[1:]:
            rv *= _sympy_to_z3_rec(var_map, child)
    elif isinstance(e, Add):
        rv = _sympy_to_z3_rec(var_map, e.args[0])

        for child in e.args[1:]:
            rv += _sympy_to_z3_rec(var_map, child)
    elif isinstance(e, Pow):
        term = _sympy_to_z3_rec(var_map, e.args[0])
        exponent = _sympy_to_z3_rec(var_map, e.args[1])

        if exponent == 0.5:
            # sqrt
            rv = z3.Sqrt(term)
        else:
            rv = term**exponent

    if rv == None:
        raise RuntimeError("Type '" + str(type(e)) + "' is not yet implemented for convertion to a z3 expresion. " + \
                            "Subexpression was '" + str(e) + "'.")

    return rv

def all_smt(s, initial_terms):
    def block_term(s, m, t):
        s.add(t != m.eval(t, model_completion=True))
    def fix_term(s, m, t):
        s.add(t == m.eval(t, model_completion=True))
    def all_smt_rec(terms):
        if z3.sat == s.check():
           m = s.model()
           yield m
           for i in range(len(terms)):
               s.push()
               block_term(s, m, terms[i])
               for j in range(i):
                   fix_term(s, m, terms[j])
               yield from all_smt_rec(terms[i:])
               s.pop()
        elif z3.unknown == s.check():
            print(z3.unknown)
    yield from all_smt_rec(list(initial_terms))


# In[300]:


def two_facet_intersection(f1_set, f2_set):
    intersection = f1_set.intersection(f2_set)
    f1_set_diff = list(f1_set - intersection)
    f1_set_diff.sort()
    f2_set_diff = list(f2_set - intersection)
    f2_set_diff.sort()
    intersection = list(intersection)
    intersection.sort()
    return intersection, f1_set_diff, f2_set_diff

def is_fangiving(K, M, symbols = [], additional_assumption = []):
    if symbols == []:
        return is_fangiving_without_indet(K, M)
    x=list(sp.symbols('x0:%d'%K.n))
    z3_vars, z3_var_map = sympy_to_z3_real_var(symbols+x)
    facet_sets = [set(facet) for facet in K.cpx]
    for i in range(len(K.cpx)):
        facet1 = facet_sets[i]
        for j in range(i+1, len(K.cpx)):
            facet2 = facet_sets[j]
            int_f = two_facet_intersection(facet1, facet2)
            N = M[:, list(np.array(int_f[0]+int_f[1]+int_f[2])-1)].nullspace()
            for k in range(len(N)):
                N[k] = sp.simplify(N[k])
            if sp.Matrix([N])[-len(int_f[2]):, :] != sp.eye(len(int_f[2])):
                print('Check dependence of facets %d, %d'%(i, j))
                return 'Check dependence of facets %d, %d'%(i, j)
            N_mat = sp.Matrix([N])
            if sp.simplify(M[:, list(np.array(int_f[0]+int_f[1]+int_f[2])-1)]*N_mat) != sp.zeros(M.shape[0], N_mat.shape[1]):
                print('Check dependence of facets %d, %d'%(i, j))
                return 'Check dependence of facets %d, %d'%(i, j)
            for k in range(len(N_mat)):
                if sp.fraction(sp.together(N_mat[k]))[1] != 1:
                    print('Check dependence of facets %d, %d'%(i, j))
                    return N_mat
            v = N[0]*x[0]
            s=z3.Solver()
            for a in additional_assumption:
                s.add(a)
            for k in range(1, len(N)):
                v = v + N[k]*x[k]
            for k in range(len(int_f[0]), len(int_f[0])+len(int_f[1])):
                s.add(sympy_to_z3_eq(z3_var_map, v[k]) <= 0)
            for k in range(len(int_f[0])+len(int_f[1]), len(int_f[0])+2*len(int_f[1])):
                s.add(sympy_to_z3_eq(z3_var_map, v[k]) >= 0)
            not_all_zero = []
            for k in range(len(N)):
                not_all_zero.append(z3_var_map[x[k].name] != 0)
            s.add(z3.Or(not_all_zero))
            check = s.check()
            if check == z3.sat:
                return (K.cpx[i], K.cpx[j], M[:, list(np.array(K.cpx[i])-1)], M[:, list(np.array(K.cpx[j])-1)], v, s.model(), N, s)
            elif check == z3.unknown:
                return (i, j, check)
    return 1


# In[290]:


class Characteristic_Pair(): # [1, ..., n] in K, M = [I_n A]
    def __init__(self, K, M):
        self.cpx = K
        self.char = sp.Matrix(M)
        if M.shape==(self.cpx.n, self.cpx.m):
            pass
        elif M.shape==(self.cpx.m, self.cpx.n):
            self.char = self.char.T
        else:
            raise Exception('The matrix size is not appropriate.')
        self._subs_stat = 1
        self.index_cpx = (np.array(self.cpx.cpx) - 1).tolist()
        self.orientation = Orientation(self.cpx)
        self._determinant = []
        self._det_ori = set()
        self._is_fangiving = None
        self.assumption = []
        D = set()
        for i in range(len(self.cpx.cpx)):
            det = sp.simplify(sp.expand(sp.det(self.char[:, self.index_cpx[i]])))
            self._determinant.append(det)
            if (det, self.orientation[i]) not in D:
                if (-det, -self.orientation[i]) not in D:
                    D.add((det, self.orientation[i]))
                    self._det_ori.add(det - self.orientation[i])          
        self._indet = set()
        for i in M:
            self._indet.update(sp.sympify(i).atoms(sp.Symbol))
        self._z3_var, self._z3_var_map = sympy_to_z3_var(self._indet)
            
    @property
    def indeterminates(self):
        return self._indet       

    @property
    def determinants(self):
        return self._determinant
        
    @property
    def is_positive(self):
        if self._det_ori == {0}:
            return 1
        else:
            return 0

    def subs(self, symbol, replaced):
        if symbol in self._indet:
            self._indet.remove(symbol)
            self._indet.update(sp.sympify(replaced).atoms(sp.Symbol))
            self.char = self.char.subs(symbol, replaced)
            D = set()
            self._det_ori = set()
            for i in range(len(self._determinant)):
                self._determinant[i] = sp.simplify(sp.expand(self._determinant[i].subs(symbol, replaced)))
                if (self._determinant[i], self.orientation[i]) not in D:
                    if (-self._determinant[i], -self.orientation[i]) not in D:
                        D.add((self._determinant[i], self.orientation[i]))
                        self._det_ori.add(self._determinant[i] - self.orientation[i])
            self._z3_var, self._z3_var_map = sympy_to_z3_var(self._indet)
            self._subs_stat = 1
            self._is_fangiving = None
            #return self.char

    def positive_map_solution(self, restart=False):
        if restart==True:
            self._subs_stat = 1
        if self._subs_stat == 1:
            self._solver = z3.Solver()
            for e in self._det_ori:
                self._solver.add(sympy_to_z3_eq(self._z3_var_map, e) == 0)
            self._sol_gen = all_smt(self._solver, self._z3_var)
            self._subs_stat = 0
        try:
            sol = next(self._sol_gen)
            new_char = self.char[:, :]
            for var in self._indet:
                new_char=new_char.subs(var, sol[self._z3_var_map[var.name]].as_long())
            return new_char
        except:
            #print('no more solution')
            return 0
            
    @property
    def is_fangiving(self):
        if self.is_positive == 0:
            print('not positive')
        if self._is_fangiving is None:
            self._is_fangiving = is_fangiving(self.cpx, self.char, list(self._indet), self.assumption)
            return self._is_fangiving
        else:
            return self._is_fangiving
    def vector(self, vertices):
        return self.char[:, [i-1 for i in vertices]]
        
    def vector_sum(self, vertices, coefficient):
        sum = sp.zeros(*self.char[:, 0].shape)
        for i in range(len(vertices)):
            sum = sum+self.char[:, vertices[i]-1]*coefficient[i]
        return sum

    def proj(self, face):
        indices1 = (np.array(face)-1).tolist()
        indices2 = (np.array(Link(self.cpx, face).vert)-1).tolist()
        relabel_M = self.char[:, indices1+indices2].rref()[0]
        return relabel_M[len(face):, len(face):]
        


# In[ ]:

def Circuits(M): # M: sympy Matrix, vector configuration
    m = len(M[0, :])
    V = set(range(m))
    d = len(M[:, 0])
    Indep = [(i, ) for i in range(m)]
    Dep = []
    Circ = []
    for i in range(2, d+2):
        Dep_next = set()
        Candi = set()
        for indep in Indep:
            Candi.update([tuple(sorted(indep+(j, ))) for j in (V-set(indep))])
        for dep in Dep:
            D = [tuple(sorted(dep+(j,))) for j in (V-set(dep))]
            Candi.difference_update(D)
            Dep_next.update(D)

        Dep = list(Dep_next)
        Indep = []
        for candi in Candi:
            r = M[:, candi].rank()
            if r == i:
                Indep.append(candi)
            else:
                Dep.append(candi)
                Circ.append(candi)
    return Circ

def Signed_Circuits(M):
    Circ = Circuits(M)
    m = len(M[0, :])
    Signed_Circ = []
    V = set(range(m))
    for circ in Circ:
        dep = M[:, circ].nullspace()[0] # minimal dependence is unique for each circuit up to scalar multiple
        pos = 0
        neg = 0
        for i in range(len(circ)):
            if dep[i] > 0:
                pos = pos + 2**(circ[i])
            elif dep[i] < 0:
                neg = neg + 2**(circ[i])
        Signed_Circ.append([pos, neg])
    return Signed_Circ

def is_fangiving_without_indet(K, M):
    facet_sets = [set(facet) for facet in K.cpx]
    SC = Signed_Circuits(M)
    for i in range(len(K.cpx)):
        facet1 = facet_sets[i]
        for j in range(i+1, len(K.cpx)):
            facet2 = facet_sets[j]
            int_f = two_facet_intersection(facet1, facet2)
            distinguished_set1 = 0
            for k in int_f[1]:
                distinguished_set1=distinguished_set1+2**(k-1)
            distinguished_set2 = 0
            for k in int_f[2]:
                distinguished_set2=distinguished_set2+2**(k-1)
            for sc in SC:
                if (sc[0] | distinguished_set1 == distinguished_set1) and (sc[1] | distinguished_set2 == distinguished_set2):
                    return K.cpx[i], K.cpx[j]
                elif (sc[1] | distinguished_set1 == distinguished_set1) and (sc[0] | distinguished_set2 == distinguished_set2):
                    return K.cpx[i], K.cpx[j]
            else:
                continue
    return 1

def initial_char(n):
    return sp.Matrix([sp.eye(n), sp.symbols('a0:%d'%n), sp.symbols('b0:%d'%n), sp.symbols('c0:%d'%n), sp.symbols('d0:%d'%n)]).T
