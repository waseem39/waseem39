
"""
Created on Mon Jul 15 15:37:37 2024

@author: Wasim
"""



from sympy import *

import cmfff as ff
# This code is usede to find the function for 

x,y,z, w, x1, y1, z1, w1, u1 =symbols('x y z w x1 y1 z1 w1 u1')

a, b, c, d, e  = symbols('a b c d e ' )

# Initial system

F = Matrix([  d * y -b*x*(z+w),
    a*(x+z)- b*y*(z+w) -(c+d)*x,
    b*x*(z+w) + d*w -e*z,
  b*y*(z+w) - (c+e+d)*w] )

Var=Matrix([x,y,z, w])

Var1= Matrix([x1,y1,z1, w1])

# Jocobian of system

A=ff.Joc(F,Var,4)


B= A.applyfunc(simplify)


B1 = B*Var

B2 = B1.applyfunc(simplify)

F1 = Matrix([b*x*w +b*x*z, b*w*y+b*y*z, -b*w*x -b*x*z, -b*w*y -b*y*z ])

F2 = (B1+ F1)

F3 = F2.applyfunc(simplify)