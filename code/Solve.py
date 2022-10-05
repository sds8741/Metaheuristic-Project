
# %%
import numpy as np
from math import *


# %%
def IDMAT(NFIX, NNOD, NDN):
    IDND = np.zeros([NDN, NNOD])
    a, b = 0, 0
    for i in range(NNOD):
        for j in range(NDN):
            if NFIX[j, i] == 0:
                a += 1
                IDND[j, i] = a
            elif NFIX[j, i] < 0:
                b -= 1
                IDND[j, i] = b
            else:
                IDND[j, i] = IDND[j, NFIX[j, i]-1]
    NEQ = np.max(IDND)
    return IDND, NEQ


# %%
def MEMDOF(IDBC, IDND, NDE, NBC, NDN):
    LM = np.zeros([NDE, NBC])
    for j in range(NBC):
        for i in range(NDE):
            a = ceil((i+1)/NDN)
            node = IDBC[a-1, j]
            k = (i+1) % NDN
            if k == 0:
                k = NDN
            LM[i, j] = IDND[k-1, node-1]
    return LM


# %%
def SEMIBAND(LM, NDE, NBC):
    A = LM.copy()
    max_LM = A.max(0)
    for i in range(NBC):
        for j in range(NDE):
            if LM[j, i] < 0:
                A[j, i] = max_LM[i]
    NSBAND = (max_LM - A.min(0) + 1).max(0)
    return NSBAND


# %%
def ELKE(NDE, IDBC, PROP, SECT, IB, RL):
    E = PROP[0, IDBC[2, IB]-1]
    NU = PROP[1, IDBC[2, IB]-1]
    A = SECT[0, IDBC[3, IB]-1]
    Iz = SECT[1, IDBC[3, IB]-1]
    Iy = SECT[2, IDBC[3, IB]-1]
    J = SECT[3, IDBC[3, IB]-1]

    EE = E*Iz/RL*np.array([[12/RL**2, 6/RL, -12/RL**2, 6/RL],
                           [6/RL, 4, -6/RL, 2],
                           [-12/RL**2, -6/RL, 12/RL**2, -6/RL],
                           [6/RL, 2, -6/RL, 4]])
    return EE


# %%
def ROTATION(COOR, IDBC, MN, NCO, NDE):
    CO = COOR[:NCO, IDBC[:2, MN]-1].T
    RL = np.sqrt(np.sum((CO[1, :] - CO[0, :])**2))
    ROT = np.eye(2)
    T = np.zeros([int(NDE), int(NDE)])
    M = 2
    for i in range(int(NDE/M)):
        dof = np.arange(M) + i*M
        s = dof[0]
        d = dof[-1]
        T[dof, s:d+1] = ROT
    return T, RL


# %%
def LOAD(EXLD, IDND, NDN, NNOD, NEQ):
    GLOAD = np.zeros([int(NEQ), 1])
    for j in range(NNOD):
        for i in range(NDN):
            if IDND[i, j] > 0:
                ind = int(IDND[i, j]) - 1
                GLOAD[ind] = GLOAD[ind] + EXLD[i, j]

    return GLOAD


# %%
def FORMKP(COOR, IDBC, PROP, SECT, LM, FEF, GLOAD, NNOD, NBC,
           NMAT, NSEC, NCO, NDN, NDE, NNE, NEQ, IFORCE):
    GLK = np.zeros([int(NEQ), int(NEQ)])

    for IB in range(NBC):
        T, RL = ROTATION(COOR, IDBC, IB, NCO, NDE)

        EE = ELKE(NDE, IDBC, PROP, SECT, IB, RL)
        LDOF = np.where(LM[:, IB] > 0)
        GDOF = LM[LDOF, IB].astype('int') - 1
        ELK = np.dot(np.dot((T.T), EE), T)

        GDOF = np.array(GDOF).reshape(-1)
        GLK_s = GLK[GDOF]
        ELK_s = ELK[LDOF]
        GLK_s[:, GDOF] = GLK_s[:, GDOF] + ELK_s[:, LDOF].squeeze()
        if IFORCE != 1:
            EFEQ = np.dot(-T.T, FEF[:, IB])
            GLOAD[GDOF] = GLOAD[GDOF] + EFEQ[LDOF].reshape(GLOAD[GDOF].shape)
        GLK[GDOF] = GLK_s
    return GLK, GLOAD


# %%
def SOLVE(GLK, GLOAD):
    return np.dot(np.linalg.inv(GLK), GLOAD)

# %%


def Execute(NNOD, NBC, NMAT, NSEC,  NNE, IFORCE, COOR, NFIX, EXLD, IDBC, PROP, SECT, FEF):

    IPR = np.array([[1, 2, 2, 2, 3, 3], [2, 2, 3, 3, 3, 6]])
    NCO = IPR[0, 0]
    NDN = IPR[1, 0]
    NDE = NDN*NNE

    IDND, NEQ = IDMAT(NFIX, NNOD, NDN)
    LM = MEMDOF(IDBC, IDND, NDE, NBC, NDN)
    GLOAD = LOAD(EXLD, IDND, NDN, NNOD, NEQ)
    GLK, GLOAD = FORMKP(COOR, IDBC, PROP, SECT, LM, FEF, GLOAD,
                        NNOD, NBC, NMAT, NSEC, NCO, NDN, NDE, NNE, NEQ, IFORCE)
    DELTA = SOLVE(GLK, GLOAD)
    return DELTA
