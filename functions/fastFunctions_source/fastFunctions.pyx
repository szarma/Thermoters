#cython: boundscheck=False, wraparound=False, nonecheck=False

import   numpy as np
cimport  numpy as np
cimport  cython

from libc.math cimport log2, log10, exp

def bindingEnergies(
               np.ndarray[np.float64_t, ndim=2] matrix,
               np.ndarray[np.int8_t  , ndim=2] seqs
               ):
    """usage: bindingEnergies(matrix, seqs)
    where matrix.shape is (L,4) and seqs.shape is (n, L),
    with n being the number of sequences"""
    
#     cdef double[::1] energies
    assert matrix.shape[0] == seqs.shape[1]
    cdef np.ndarray[np.float64_t, ndim=1] energies
    cdef unsigned int ns, L, i,j                    
    
    ns = len(seqs)
    L  = len(seqs[0])
    energies = np.zeros(ns)
    
    for i in range(ns):
        for j in range(L):
            energies[i] += matrix[j,seqs[i,j]]
            
    return energies

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def getDiNu(long p1,
            long b1, 
            long p2, 
            long b2, 
            long n1, 
            long minSpacer, 
            long n2, 
            np.ndarray[np.int8_t  , ndim=2] sequences, 
            long nSpacer, 
            ):

    cdef int minMsize#,i,iS,x0,js    
    cdef Py_ssize_t nSeq = sequences.shape[0]
    cdef Py_ssize_t seqL = sequences.shape[1]
    cdef Py_ssize_t js,i,x0,iS
    
    minMsize = minSpacer+n1+n2

    diNu = np.zeros((nSpacer,seqL-minMsize-nSpacer+1,nSeq),dtype=np.intc)#.astype(np.int32)
    cdef int[:,:,::1] diNu_view = diNu
    
    if p1>=n1: p1 -= nSpacer//2
    if p2>=n1: p2 -= nSpacer//2

    for iS in xrange(nSpacer):
        x0 = nSpacer-iS
        for i in xrange(x0,seqL-minMsize-iS+1):
            for js in xrange(nSeq):
                if sequences[js][i+p1]==b1 and sequences[js][i+p2]==b2:
                    diNu_view[iS,i-x0,js] = 1
#                     diNu[iS,i-x0,js] = 1
        if p1>=n1: p1 += 1
        if p2>=n1: p2 += 1

    return diNu

from numpy import zeros
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def tensum(
    double[::1] a_,
    int[:,:,:,::1] B_):
    cdef Py_ssize_t nInter_ = B_.shape[0]
    cdef Py_ssize_t nSeq_   = B_.shape[1]
    cdef Py_ssize_t nPos_   = B_.shape[2]
    cdef Py_ssize_t nSpacer_= B_.shape[3]
    
    out = zeros((nSeq_,nPos_,nSpacer_),dtype=np.double)
    cdef double[:,:,::1] out_view = out#zeros((nSeq_,nPos_,nSpacer_))
    cdef int i,j,k,l
    for i in xrange(nSeq_):
        for j in xrange(nPos_):
            for k in xrange(nSpacer_):
                for l in xrange(nInter_):
                    out_view[i,j,k] += a_[l]*B_[l,i,j,k]
#                     out[i,j,k] += a_[l]*B_[l,i,j,k]
    
    return out

