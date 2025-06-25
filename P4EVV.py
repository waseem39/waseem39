# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 01:05:28 2024

@author: Wasim
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 21:10:01 2024

@author: Wasim
"""

from sympy import *

import cmfff as ff

x, y, z, w = symbols('x y z w')

a, b, c , d , e ,t = symbols('a, b, c, d, e ,t')
 
Var = (x,y,z,w)
F1 =Matrix([0,(a*c)/(b*(c+a)),c*t/a,0])

F = Matrix([  d * y -b*x*(z+w),
    a*(x+z)- b*y*(z+w) -(c+d)*x,
    
    b*x*(z+w) + d*w -e*z,
  b*y*(z+w) - (c+e+d)*w] )

 
A = ff.Joc(F, Var, 4)

B = A.subs({x:F1[0], y:F1[1], z: F1[2], w: F1[3]})

B1 = B.subs({d:0,e:0})

ev =B1.eigenvals()

# evv =B1.eigenvects()


# BT = B1.transpose()

# evv2 = BT.eigenvects()




