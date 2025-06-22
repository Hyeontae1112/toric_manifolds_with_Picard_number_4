#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
from itertools import combinations, product, permutations
from tqdm import tqdm
import json
def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data
def text(result, path):
    t = open(path, mode='a', encoding='utf-8')
    for sol in result:
        t.write(str(sol) + '\n')
    t.close()

# In[43]:


def Cpx(K2):
    K=[]
    for k in K2:
        K.append([])
        i=1
        while 2**(i-1)<k+1:
            if 2**(i-1)|k == k:
                K[-1].append(i)
            i+=1
    return K

def Cpx_bin(K):
    K2=[]
    for k in K:
        z=0
        for i in k:
            z=z+2**(i-1)
        K2.append(z)
    return K2

def Vert(K2):
    a=K2[0]
    for i in K2:
        a=a|i
    return (Cpx([a])[0], a)

def Face_Poset(K, n):
    Faces2=[]
    f_vector = [0]*n
    for d in range(1, n+1):
        d_faces=set()
        for facet_idx in range(len(K)):
            C=combinations(K[facet_idx], d)
            for c in C:
                d_faces.add(Cpx_bin([c])[0])
        d_faces=list(d_faces)
        f_vector[d-1] = len(d_faces)
        Faces2.extend(d_faces)
    return Faces2, f_vector

def Z2n(n):
    V=[]
    N=range(n)
    for i in range(1,n+1):
        C=combinations(N,i)
        for c in C:
            v=[]
            for j in range(n):
                if j in c:
                    v.append(1)
                    continue
                v.append(0)
            V.append(v)
    return V

def Z2n_bin(n):
    V=[]
    E=[2**(i) for i in range(n-1, -1, -1)]
    N=range(n)
    for i in range(1,n+1):
        C=combinations(E,i)
        for c in C:
            V.append(sum(c))
    return V    

def relabel(K, V, W): # K=cpx, V->W
    L=[]
    for facet in K:
        L.append([])
        for v in facet:
            L[-1].append(W[V.index(v)])
        L[-1].sort()
    L.sort()
    return L

def Minimal_non_faces(cpx_bin, m, n, vert):
    MNF = []
    mnf_vector = [0]*n
    for d in range(2, n+2):
        C = combinations(vert, d)
        for c in C:
            nf = Cpx_bin([c])[0]
            for mnf in MNF:
                if mnf|nf == nf:
                    break
            else:
                for f in cpx_bin:
                    if nf|f == f:
                        break
                else:
                    MNF.append(nf)
                    mnf_vector[d-2] += 1
    return MNF, mnf_vector


# In[13]:


class Simplicial_Complex():
    def __init__(self, K):
        if type(K[0]) == int:
            self.cpx = Cpx(K)
            self.cpx.sort()
            self.cpx_bin = Cpx_bin(self.cpx)
        else:
            self.cpx = K
            for i in self.cpx:
                i.sort()
            self.cpx.sort()
            self.cpx_bin = Cpx_bin(K)
            
        self.n = len(self.cpx[0])
        self.dim = self.n-1
        self.vert, self.vert_bin = Vert(self.cpx_bin)
        self.m = len(self.vert)
        self.p=self.m-self.n
        self.min_non_faces = None
        self._faces=None
        self._faces_bin=None
        self._cofacets=None
        self._cofacets_bin=None
        self._mnf=None
        self._mnf_bin=None
        self._f_vector=None
        self._mnf_vector=None
        self._symm=None
        
    def Cofacets(self):
        self._cofacets_bin=[]
        for i in self.cpx_bin:
            self._cofacets_bin.append(i^self.vert_bin)
        self._cofacets = Cpx(self._cofacets_bin)        
    @property
    def cofacets_bin(self):
        if self._cofacets_bin is None:
            self.Cofacets()
        return self._cofacets_bin
    @property
    def cofacets(self):
        if self._cofacets is None:
            self.Cofacets()
        return self._cofacets
    
    def Faces(self):
        self._faces_bin, self._f_vector = Face_Poset(self.cpx, self.n)
        self._faces=Cpx(self._faces_bin)
    @property
    def faces_bin(self):
        if self._faces_bin is None:
            self.Faces()
        return self._faces_bin
    @property
    def faces(self):
        if self._faces is None:
            self.Faces()
        return self._faces
    @property
    def f_vector(self):
        if self._f_vector is None:
            self.Faces()
        return self._f_vector
    
    def is_seed(self):
        C=combinations(self.vert, 2)
        for c in C:
            for i in self.cpx_bin:
                if (2**(c[0]-1))|i != i and (2**(c[1]-1))|i != i:
                    break
            else:
                for i in self.cpx_bin:
                    if (2**(c[0]-1))|i == i and (2**(c[1]-1))|i == i:
                        return 0
        return 1
    
    def MNF(self):
        self._mnf_bin, self._mnf_vector = Minimal_non_faces(self.cpx_bin, self.m, self.n, self.vert)
        self._mnf=Cpx(self._mnf_bin)
    @property
    def mnf_bin(self):
        if self._mnf_bin is None:
            self.MNF()
        return self._mnf_bin
    @property
    def mnf(self):
        if self._mnf is None:
            self.MNF()
        return self._mnf
    @property
    def mnf_vector(self):
        if self._mnf_vector is None:
            self.MNF()
        return self._mnf_vector

    @property
    def symmetries(self):
        if self._symm is None:
            self._symm = symmetries(self)
        return self._symm
        


# In[25]:


def CM_bin(K, l):
    #vertex relabeling
    original_vert = Cpx([K.cpx_bin[0]])[0]+Cpx([K.vert_bin ^ K.cpx_bin[0]])[0]
    K = Simplicial_Complex(relabel(K.cpx, original_vert, K.vert))
    relabeled_vert=[]
    for v in K.vert:
        relabeled_vert.append(original_vert.index(v))
    
    m=K.m
    n=K.n
    if m == n:
        return [np.eye(n)]
    result = []
    Z=Z2n_bin(l)
    L=np.zeros(m, dtype=int)
    L[:n] = Z[:n]
    column_index = n
    vertex = n+1
    candidate_vectors=[0]*m
    check=1
    faces=[]
    for f in K.faces:
        if len(f)>1:
            faces.append(np.array(f)-1)
    while column_index != n-1:
        if check==1:
            candidate_vectors[column_index]=Z[:]
            for face in faces:
                if face[-1]==column_index:
                    s=np.bitwise_xor.reduce(L[face[:-1]])
                    if s in candidate_vectors[column_index]:
                            candidate_vectors[column_index].remove(s)
                            
        if candidate_vectors[column_index]==[]:
            vertex = vertex - 1
            column_index = column_index - 1
            check=0
            continue
        L[column_index] = candidate_vectors[column_index][0]
        del candidate_vectors[column_index][0]
        if vertex == m:
            check = 0
            result.append(L[relabeled_vert])
            continue
        vertex = vertex + 1
        column_index = column_index + 1
        check = 1
    return result


# In[87]:


def CM(K, l):
    cm_bin = CM_bin(K,l)
    result = [[[] for j in range(K.m)] for i in cm_bin]
    for i in range(len(cm_bin)):
        result[i] = ((cm_bin[i].reshape(-1,1) & (2**np.arange(K.n-1, -1, -1))) != 0).astype(int)
    return result


# In[29]:


# def CM(K, l):
#     #vertex relabeling
#     original_vert = Cpx([K.cpx_bin[0]])[0]+Cpx([K.vert_bin ^ K.cpx_bin[0]])[0]
#     K = Simplicial_Complex(relabel(K.cpx, original_vert, K.vert))
#     relabeled_vert=[]
#     for v in K.vert:
#         relabeled_vert.append(original_vert.index(v))
    
#     m=K.m
#     n=K.n
#     if m == n:
#         return [np.eye(n)]
#     result = []
#     Z=Z2n(l)
#     L=np.zeros((m,l), dtype=int)
#     L[:n] = Z[:n]
#     column_index = n
#     vertex = n+1
#     candidate_vectors=[0]*m
#     check=1
#     faces=[]
#     for f in K.faces:
#         faces.append(np.array(f)-1)
#     while column_index != n-1:
#         if check==1:
#             candidate_vectors[column_index]=Z[:]
#             for face in faces:
#                 if face[-1]==column_index and len(face) > 1:
#                     s=sum(L[face[:-1]])%2
#                     for c in range(len(candidate_vectors[column_index])):
#                         if (candidate_vectors[column_index][c] == s).all():
#                             del candidate_vectors[column_index][c]
#                             break
                            
#         if candidate_vectors[column_index]==[]:
#             vertex = vertex - 1
#             column_index = column_index - 1
#             check=0
#             continue
#         L[column_index] = candidate_vectors[column_index][0]
#         del candidate_vectors[column_index][0]
#         if vertex == m:
#             check = 0
#             result.append(L.copy()[relabeled_vert])
#             continue
#         vertex = vertex + 1
#         column_index = column_index + 1
#         check = 1
#     return result


# In[31]:


def is_colorable(K, l): # if 1 is returned, then the Buchstaber number of K is larger than or equal to m-l
    #vertex relabeling
    original_vert = Cpx([K.cpx_bin[0]])[0]+Cpx([K.vert_bin ^ K.cpx_bin[0]])[0]
    K = Simplicial_Complex(relabel(K.cpx, original_vert, K.vert))
    
    m=K.m
    n=K.n
    Z=Z2n_bin(l)
    L=np.zeros(m, dtype=int)
    L[:n] = Z[:n]
    column_index = n
    vertex = n+1
    candidate_vectors=[0]*m
    check=1
    faces=[]
    for f in K.faces:
        if len(f)>1:
            faces.append(np.array(f)-1)
    while column_index != n-1:
        if check==1:
            candidate_vectors[column_index]=Z[:]
            for face in faces:
                if face[-1]==column_index:
                    s=np.bitwise_xor.reduce(L[face[:-1]])
                    if s in candidate_vectors[column_index]:
                            candidate_vectors[column_index].remove(s)

                            
        if candidate_vectors[column_index]==[]:
            vertex = vertex - 1
            column_index = column_index - 1
            check=0
            continue
        L[column_index] = candidate_vectors[column_index][0]
        del candidate_vectors[column_index][0]
        if vertex == m:
            check = 0
            return 1
        vertex = vertex + 1
        column_index = column_index + 1
        check = 1
    return 0


# In[6]:


def Real_Buchstaber_number(K):
    for l in range(K.n, K.m):
        if is_colorable(K, l)==1:
            return K.m-l
    return 0


# In[7]:


def Wedge(K, v):
    
    Wedge_K = []
    v_bin = 2**(v-1)
    new_vert_bin = 2**v
    
    if v_bin | K.vert_bin != K.vert_bin:
        return 0
    
    for facet_bin in K.cpx_bin:
        f1 = facet_bin >> v << v
        f2 = facet_bin ^ f1
        new_facet_bin = (f1 << 1) + f2

        if v_bin | facet_bin == facet_bin:
            Wedge_K.append(new_facet_bin + 2*v_bin)
            
        else:
            Wedge_K.append(new_facet_bin + v_bin)
            Wedge_K.append(new_facet_bin + 2*v_bin)
            
    return Simplicial_Complex(Wedge_K)


# In[8]:


def J_construction(K, J):
    Wedge_K = K.cpx_bin
    v=1
    for j in J:
        K = Wedge_K
        Wedge_K = []
        v_bin = 2**(v-1)
        for facet_bin in K:
            f1 = facet_bin >> v << v
            f2 = facet_bin ^ f1
            new_facet_bin = (f1 << (j-1)) + f2
            if v_bin | facet_bin == facet_bin:
                for i in range(2, j+1):
                    new_facet_bin = new_facet_bin + (v_bin << (i-1))
                Wedge_K.append(new_facet_bin)
            else:
                for i in range(1, j+1):
                    new_facet_bin = new_facet_bin + (v_bin << (i-1))
                for i in range(1, j+1):
                    Wedge_K.append(new_facet_bin - (v_bin << (i-1)))
        v += j    
    return Simplicial_Complex(Wedge_K)


# In[9]:


def Suspension(K):
    
    Susp_K = []
    new_vert_bin = 2**(max(K.vert) + 1)
    
    for facet_bin in K.cpx_bin:
        f = facet_bin << 1
        Susp_K.append(1 + f)
        Susp_K.append(f + new_vert_bin)
        
    return Simplicial_Complex(Susp_K)


# In[99]:


def DCM_bin(K):
    #vertex relabeling
    original_vert = Cpx([K.cpx_bin[0]])[0]+Cpx([K.vert_bin ^ K.cpx_bin[0]])[0]
    K = Simplicial_Complex(relabel(K.cpx, original_vert, K.vert))
    relabeled_vert=[]
    for v in K.vert:
        relabeled_vert.append(original_vert.index(v))
    
    m=K.m
    n=K.n
    p=m-n
    result = []
    Z=Z2n_bin(p)
    L=np.zeros(m, dtype=int)
    L[n:] = Z[:p]
    row_index = n-1
    vertex = n
    candidate_vectors=[0]*m
    check=1
    D_facets = K.cofacets
    D_K = Simplicial_Complex(D_facets)
    D_faces=[]
    for D_f in D_K.faces:
        if len(D_f) > 1:
            D_faces.append(np.array(D_f)-1)
    while row_index != n:
        if check==1:
            candidate_vectors[row_index]=Z[:]
            for D_face in D_faces:
                if D_face[0]==row_index and len(D_face) > 1:
                    s=np.bitwise_xor.reduce(L[D_face[1:]])
                    if s in candidate_vectors[row_index]:
                        candidate_vectors[row_index].remove(s)
                            
        if candidate_vectors[row_index]==[]:
            vertex = vertex + 1
            row_index = row_index + 1
            check=0
            continue
        L[row_index] = candidate_vectors[row_index][0]
        del candidate_vectors[row_index][0]
        if vertex == 1:
            check = 0
            result.append(L[relabeled_vert])
            continue
        vertex = vertex - 1
        row_index = row_index - 1
        check = 1
    return result


# In[102]:


def DCM(K):
    dcm_bin = DCM_bin(K)
    result = [0]*len(dcm_bin)
    for i in range(len(dcm_bin)):
        result[i] = ((dcm_bin[i].reshape(-1,1) & (2**np.arange(K.p-1, -1, -1))) != 0).astype(int)
    return result


# In[106]:


def IDCM_bin(K): 
    #vertex relabeling
    original_vert = Cpx([K.cpx_bin[0]])[0]+Cpx([K.vert_bin ^ K.cpx_bin[0]])[0]
    K = Simplicial_Complex(relabel(K.cpx, original_vert, K.vert))
    relabeled_vert=[]
    for v in K.vert:
        relabeled_vert.append(original_vert.index(v))
    
    m=K.m
    n=K.n
    p=m-n
    result = []
    Z=Z2n_bin(p)
    L=np.zeros(m, dtype=int)
    L[n:] = Z[:p]
    Z=Z[p:]
    row_index = n-1
    vertex = n
    candidate_vectors=[0]*m
    unused_vectors = [0]*m
    unused_vectors[n] = Z[:]
    check=1
    D_facets = K.cofacets
    D_K = Simplicial_Complex(D_facets)
    D_faces=[]
    for D_f in D_K.faces:
        if len(D_f) > 1:
            D_faces.append(np.array(D_f)-1)
    while row_index != n:
        if check==1:
            candidate_vectors[row_index]=unused_vectors[row_index + 1][:]
            for D_face in D_faces:
                if D_face[0]==row_index:
                    s=np.bitwise_xor.reduce(L[D_face[1:]])
                    if s in candidate_vectors[row_index]:
                        candidate_vectors[row_index].remove(s)
                            
        if candidate_vectors[row_index]==[]:
            vertex = vertex + 1
            row_index = row_index + 1
            check=0
            continue
        L[row_index] = candidate_vectors[row_index][0]
        unused_vectors[row_index] = unused_vectors[row_index + 1][:]
        unused_vectors[row_index].remove(candidate_vectors[row_index][0])
        del candidate_vectors[row_index][0]
        if vertex == 1:
            check = 0
            result.append(L[relabeled_vert])
            continue
        vertex = vertex - 1
        row_index = row_index - 1
        check = 1
    return result


# In[107]:


def IDCM(K):
    idcm_bin = IDCM_bin(K)
    result = [0]*len(idcm_bin)
    for i in range(len(idcm_bin)):
        result[i] = ((idcm_bin[i].reshape(-1,1) & (2**np.arange(K.p-1, -1, -1))) != 0).astype(int)
    return result


# In[37]:


def is_IDCM(K): 
    #vertex relabeling
    original_vert = Cpx([K.cpx_bin[0]])[0]+Cpx([K.vert_bin ^ K.cpx_bin[0]])[0]
    K = Simplicial_Complex(relabel(K.cpx, original_vert, K.vert))
    relabeled_vert=[]
    for v in K.vert:
        relabeled_vert.append(original_vert.index(v))
    
    m=K.m
    n=K.n
    p=m-n
    result = []
    Z=Z2n_bin(p)
    L=np.zeros(m, dtype=int)
    L[n:] = Z[:p]
    Z=Z[p:]
    row_index = n-1
    vertex = n
    candidate_vectors=[0]*m
    unused_vectors = [0]*m
    unused_vectors[n] = Z[:]
    check=1
    D_facets = K.cofacets
    D_K = Simplicial_Complex(D_facets)
    D_faces=[]
    for D_f in D_K.faces:
        if len(D_f) > 1:
            D_faces.append(np.array(D_f)-1)
    while row_index != n:
        if check==1:
            candidate_vectors[row_index]=unused_vectors[row_index + 1][:]
            for D_face in D_faces:
                if D_face[0]==row_index:
                    s=np.bitwise_xor.reduce(L[D_face[1:]])
                    if s in candidate_vectors[row_index]:
                        candidate_vectors[row_index].remove(s)
                            
        if candidate_vectors[row_index]==[]:
            vertex = vertex + 1
            row_index = row_index + 1
            check=0
            continue
        L[row_index] = candidate_vectors[row_index][0]
        unused_vectors[row_index] = unused_vectors[row_index + 1][:]
        unused_vectors[row_index].remove(candidate_vectors[row_index][0])
        del candidate_vectors[row_index][0]
        if vertex == 1:
            check = 0
            return 1
        vertex = vertex - 1
        row_index = row_index - 1
        check = 1
    return 0


# In[229]:


def is_iso(K1, K2):
    if K1.f_vector != K2.f_vector:
        return 0
    if K1.mnf_vector != K2.mnf_vector:
        return 0
    color1 = [[] for i in range(K1.m)]
    color2 = [[] for i in range(K2.m)]
    mnf1 = K1.mnf
    mnf2 = K2.mnf

    for v in K1.vert:
        for m in K1.mnf_bin:
            if 2**(v-1) | m == m:
                color1[v-1].append(m.bit_count())
        color1[v-1].sort()
    for v in K2.vert:
        for m in K2.mnf_bin:
            if 2**(v-1) | m == m:
                color2[v-1].append(m.bit_count())
        color2[v-1].sort()
    
    partition1 = []
    partition2 = []
    for i in range(len(color1)):
        if color1[i] not in partition1:
            partition1.append(color1[i])
            partition2.append([])
    for i in range(len(partition1)):
        c = partition1[i]
        partition1[i] = []
        for j in range(len(color1)):
            if color1[j] == c:
                partition1[i].append(j+1)
        for j in range(len(color2)):
            if color2[j] == c:
                partition2[i].append(j+1)
                
    color1.sort()
    color2.sort()
    if color1 != color2:
        return 0
    K2_mnf = []
    start = 0
    end = 0
    for f in K2.mnf_vector:
        start = end
        end = start+f
        K2_mnf.append(K2.mnf[start:end])
        K2_mnf[-1].sort() 
    R = []
    for i in range(len(partition1)):
        R.append(permutations(partition1[i]))
    P = product(*R)
    W = [v for sublist in partition2 for v in sublist]
    for p in P:
        V = [v for sublist in p for v in sublist]
        relabeled_K1_mnf = []
        start=0
        end=0
        for f in K1.mnf_vector:
            start = end
            end = start+f
            relabeled_K1_mnf.append(relabel(K1.mnf[start:end],V,W))
        if relabeled_K1_mnf == K2_mnf:
            return 1
    return 0


# In[1]:


def Orientation(K):
    O = [0]*len(K.cpx)
    ridge_index = []
    assigned_facet_index = set()
    unassigned_facet_index = set(range(len(K.cpx)))
    single_ridge_index = set()
    
    for i in range(len(K.faces)):
        if len(K.faces[i]) == K.n - 1:
            ridge_index.append(i)
    
    B = [set() for i in K.cpx]
    for i in ridge_index:
        for j in range(len(K.cpx)):
            if K.faces_bin[i] | K.cpx_bin[j] == K.cpx_bin[j]:
                B[j].add(i)
    
    O[0] = 1
    unassigned_facet_index.remove(0)
    assigned_facet_index.add(0)
    single_ridge_index.update(B[0])
    
    while len(unassigned_facet_index) != 0:
        for i in single_ridge_index:
            for j in unassigned_facet_index:
                if K.faces_bin[i] | K.cpx_bin[j] == K.cpx_bin[j]:
                    f_index1 = j
                    break
            for j in assigned_facet_index:
                if K.faces_bin[i] | K.cpx_bin[j] == K.cpx_bin[j]:
                    f_index2 = j
                    break
            for v in K.cpx[f_index1]:
                if v not in K.faces[i]:
                    v1 = v
                    break
            for v in K.cpx[f_index2]:
                if v not in K.faces[i]:
                    v2 = v
                    break
            O[f_index1] = O[f_index2] * (-1)**(K.cpx[f_index1].index(v1) + K.cpx[f_index2].index(v2) + 1)                    
            for r in B[f_index1]:
                if r in single_ridge_index:
                    single_ridge_index.remove(r)
                else:
                    single_ridge_index.add(r)
            assigned_facet_index.add(f_index1)
            unassigned_facet_index.remove(f_index1)
            break
    return O

def Orientation_check(K, O):
    ridges = []
    ridges_gen = []
    for f in K.faces:
        if len(f) == K.n - 1:
            ridges.append(f)
            ridges_gen.append(0)
                
    for i in range(len(K.cpx)):
        for j in range(len(K.cpx[i])):
            R = K.cpx[i][:j] + K.cpx[i][j+1:]
            R_index = ridges.index(R)
            ridges_gen[R_index] = ridges_gen[R_index] + (-1)**j*O[i]
    return ridges_gen

def Stellar_subdivision(K, f):
    new_v_bin = 2**K.m
    f_index = K.faces.index(f)
    partial_f = []
    for i in range(0, len(f)):
        partial_f.append(f[:i] + f[i+1:])
    partial_f_bin = Cpx_bin(partial_f)
    
    Ss_K = []
    for facet_bin in K.cpx_bin:
        if K.faces_bin[f_index] | facet_bin != facet_bin:
            Ss_K.append(facet_bin)
            
    for facet_bin in K.cpx_bin:
        if K.faces_bin[f_index] | facet_bin == facet_bin:
            for p_f_bin in partial_f_bin:
                Ss_K.append(new_v_bin + facet_bin - K.faces_bin[f_index] + p_f_bin)
    
    return Simplicial_Complex(Ss_K)

def Link(K, face, relabeling=False):
    L = []
    f_bin=K.faces_bin[K.faces.index(face)]
    for facet_bin in K.cpx_bin:
        if f_bin | facet_bin == facet_bin:
            L.append(f_bin^facet_bin)
    L = Simplicial_Complex(L)
    if relabeling==True:
        return Simplicial_Complex(relabel(L.cpx, L.vert, list(range(1, L.m+1))))
    return L

def symmetries(K):
    symm = []
    color = [[] for i in range(K.m)]
    mnf = K.mnf

    for v in K.vert:
        for m in K.mnf_bin:
            if 2**(v-1) | m == m:
                color[v-1].append(m.bit_count())
        color[v-1].sort()
    
    partition = []
    for i in range(len(color)):
        if color[i] not in partition:
            partition.append(color[i])
    for i in range(len(partition)):
        c = partition[i]
        partition[i] = []
        for j in range(len(color)):
            if color[j] == c:
                partition[i].append(j+1)
    color.sort()
    
    K_mnf = []
    start = 0
    end = 0
    for f in K.mnf_vector:
        start = end
        end = start+f
        K_mnf.append(K.mnf[start:end])
        K_mnf[-1].sort() 
    R = []
    for i in range(len(partition)):
        R.append(permutations(partition[i]))
    P = product(*R)
    W = [v for sublist in partition for v in sublist]
    for p in P:
        V = [v for sublist in p for v in sublist]
        relabeled_K_mnf = []
        start=0
        end=0
        for f in K.mnf_vector:
            start = end
            end = start+f
            relabeled_K_mnf.append(relabel(K.mnf[start:end], W, V))
        if relabeled_K_mnf == K_mnf:
            symm.append([V[W.index(v)] for v in K.vert])
    return symm
