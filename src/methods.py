#!/usr/bin/env python

from mpmath import *

def p_deriv(p, order=0):
    ''' calculate mpmath polynomial derivative coefficients '''
    p_ = p.copy()
    for i in range(order):
        for j, n in enumerate(arange(p_.rows - i)[::-1]):
            p_[j] *= n
    return p_ if order==0 else p_[:-order]

def import_profile(file):
    ''' import scan lengths (um) from taper profile file '''
    with open(file, mode='r') as f: 
        x = [int(x) for x in f.read().split()]
    return [mpf(x_)/mpf('1e6') for x_ in x]

def init_propagate(B):
    ''' perform initial propagation of recursive values for fpga '''
    b = B.copy()
    for i in range(len(b)-1):
        for j in range(len(b)-(i+2)):
            b[j+1] += b[j]
    return b

def propagate(b):
    ''' perform single propagation of recursive values '''
    for i in range(len(b)-1):
        b[len(b)-(i+1)] += b[len(b)-(i+2)]
    return b

class jerk_matrices:
    def __init__(self, dt, dx, bits):
        self.dt = dt
        self.dx = dx
        self.bits = bits
    
    def initial(self, t, x):
        ''' Initial stage movement starting from standstill '''
        A = matrix([[0,       0,      0,      0,    0, 1],    #       x_0 = 0
                    [0,       0,      0,      0,    1, 0],    #   dx/dt_0 = 0
                    [0,       0,      0,      2,    0, 0],    # d2x/dt2_0 = 0
                    [t**5,    t**4,   t**3,   t**2, t, 1],    #       x_t = x
                    [5*t**4,  4*t**3, 3*t**2, 2*t,  1, 0],    #   dx/dt_t = 0
                    [60*t**2, 24*t,   6,      0,    0, 0]])   # d3x/dt3_t = 0

        return A**-1 * matrix([0, 0, 0, x, 0, 0])

    def final(self, t, d2x_0):
        ''' Final stage movement ending at standstill '''
        A = matrix([[0,       0,       0,      0,    0, 1],    #       x_0 = 0
                    [0,       0,       0,      0,    1, 0],    #   dx/dt_0 = 0
                    [0,       0,       0,      2,    0, 0],    # d2x/dt2_0 = d2x_0
                    [0,       0,       6,      0,    0, 0],    # d3x/dt3_0 = 0
                    [5*t**4,  4*t**3,  3*t**2, 2*t,  1, 0],    #   dx/dt_t = 0
                    [20*t**3, 12*t**2, 6*t,    2,    0, 0]])   # d2x/dt2_t = 0

        return A**-1 * matrix([0, 0, d2x_0, 0, 0, 0])

    def mid(self, t, d2x_0):
        ''' Stage movement during pull '''
        A = matrix([[0,       0,      0,      0,   0, 1],    #       x_0 = 0
                    [0,       0,      0,      0,   1, 0],    #   dx/dt_0 = 0
                    [0,       0,      0,      2,   0, 0],    # d2x/dt2_0 = d2x_0
                    [0,       0,      6,      0,   0, 0],    # d3x/dt3_0 = 0
                    [5*t**4,  4*t**3, 3*t**2, 2*t, 1, 0],    #   dx/dt_t = 0
                    [60*t**2, 24*t,   6,      0,   0, 0]])   # d3x/dt3_t = 0

        b = matrix([0, 0, d2x_0, 0, 0, 0])

        return A**-1 * b

    def hogan(self, t, x):
        ''' minimum jerk solution based on flash1985 '''
        return matrix([x*6/t**5, -x*15/t**4, x*10/t**3, 0, 0, 0])

    def jerk_cost(self, p, t):
        ''' integrate jerk cost function from 0 to t '''
        return 12*(60*p[0]**2*t**5 + 60*p[1]*p[0]*t**4 + (16*p[1]**2 + 20*p[2]*p[0])*t**3 + 12*p[2]*p[1]*t**2 + 3*p[2]**2*t)

    def to_recursive(self, p):
        ''' convert polynomial coefficients to fpga recursive values '''
        A_inv = matrix([[ 120,    0,   0,   0,  0,  0],
                        [-240,   24,   0,   0,  0,  0],
                        [ 150,  -36,   6,   0,  0,  0],
                        [ -10,   14,  -6,   2,  0,  0],
                        [  -9,   -1,   1,  -1,  1,  0],
                        [   0,    0,   0,   0,  0,  1]])

        b = A_inv * p
        for i, b_ in enumerate(b):
            b[i] = nint(b_*2**self.bits)
        return b

    def segmented_initial(self, t, x, dx_t):
        ''' Initial stage movement starting from standstill '''
        A = matrix([[0,       0,       0,      0,    0, 1],    #       x_0 = 0
                    [0,       0,       0,      0,    1, 0],    #   dx/dt_0 = 0
                    [0,       0,       0,      2,    0, 0],    # d2x/dt2_0 = 0
                    [t**5,    t**4,    t**3,   t**2, t, 1],    #       x_t = x
                    [5*t**4,  4*t**3,  3*t**2, 2*t,  1, 0],    #   dx/dt_t = dx_t
                    [20*t**3, 12*t**2, 6*t,    2,    0, 0]])   # d2x/dt2_t = 0

        return A**-1 * matrix([0, 0, 0, x, dx_t, 0])

    def segmented_mid(self, t, x):
        ''' Stage movement at constant velocity '''
        A = matrix([[0,       0,       0,      0,    0, 1],    #       x_0 = 0
                    [0,       0,       0,      0,    1, 0],    #   dx/dt_0 = x/t
                    [0,       0,       0,      2,    0, 0],    # d2x/dt2_0 = 0
                    [t**5,    t**4,    t**3,   t**2, t, 1],    #       x_t = x
                    [5*t**4,  4*t**3,  3*t**2, 2*t,  1, 0],    #   dx/dt_t = x/t
                    [20*t**3, 12*t**2, 6*t,    2,    0, 0]])   # d2x/dt2_t = 0

        return A**-1 * matrix([0, x/t, 0, x, x/t, 0])


    def segmented_final(self, t, x, dx_t):
        ''' Final stage movement starting from speed dx_t '''
        A = matrix([[0,       0,       0,      0,    0, 1],    #       x_0 = 0
                    [0,       0,       0,      0,    1, 0],    #   dx/dt_0 = dx_t
                    [0,       0,       0,      2,    0, 0],    # d2x/dt2_0 = 0
                    [t**5,    t**4,    t**3,   t**2, t, 1],    #       x_t = x
                    [5*t**4,  4*t**3,  3*t**2, 2*t,  1, 0],    #   dx/dt_t = 0
                    [20*t**3, 12*t**2, 6*t,    2,    0, 0]])   # d2x/dt2_t = 0

        return A**-1 * matrix([0, dx_t, 0, x, 0, 0])
