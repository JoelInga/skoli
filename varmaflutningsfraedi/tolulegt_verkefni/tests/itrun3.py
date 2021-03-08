from scipy.optimize import fsolve
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt
import numpy as np

# n   : Fjöldi eininga
# L   : Lengd svæðis
# rho : Eðlismassi
# c   : Eðlisvarmi
# k   : Varmaleiðnistuðull
# h   : Varmaburðarstuðull (tvö gildi, vinstri og hægri)
# Tinf: Hitastig í fjarlægð (tvö gildi, vinstri og hægri)

T_inf = np.array([95,18])
n = 100


def fem(n, L, rho, c, k, h, Tinf):
    dx = L/n
    M = sp.dok_matrix((n+1,n+1))
    A = sp.dok_matrix((n+1,n+1))
    b = np.zeros(n+1)
    for i in range(n):
        M[i:i+2,i:i+2] += np.array([[2,1],[1,2]]) * rho * c * dx / 6.0
        A[i:i+2,i:i+2] += np.array([[1,-1],[-1,1]]) * k / dx
    for i,j in enumerate([0,n]):
        A[j,j] += h[i]
        b[j] += h[i] * Tinf[i]
    
    return A.tocsc(), b, M.tocsc()


def femaxi(n, ri, ro, rho, c, k, h, Tinf):
    r = np.linspace(ri,ro,n+1)
    M = sp.dok_matrix((n+1,n+1))
    A = sp.dok_matrix((n+1,n+1))
    b = np.zeros(n+1)
    for i in range(n):
        M[i:i+2,i:i+2] += np.array([[3*r[i]+r[i+1],r[i]+r[i+1]],
            [r[i]+r[i+1],r[i]+3*r[i+1]]]) * rho * c * np.pi * (r[i+1]-r[i]) / 6.0
        A[i:i+2,i:i+2] += np.array([[1,-1],[-1,1]]) \
            * k * np.pi * (r[i+1]+r[i])/(r[i+1]-r[i])
    for i,j in enumerate([0,n]):
        A[j,j] += h[i] * 2*np.pi*r[j]
        b[j] += h[i] * Tinf[i] * 2*np.pi*r[j]
    
    return A.tocsc(), b, M.tocsc()