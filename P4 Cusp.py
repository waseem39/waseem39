# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:41:03 2024

@author: Wasim
"""

from sympy import*

import cmfff as ff

x, y, z, w,  x1, y1, z1 ,w1 = symbols('x y z w x1, y1, z1, w1', reral =True)

a, b, c , d , e ,t = symbols('a, b, c, d, e ,t')
 
Var = (x,y,z,w)

Var1 =  (x1,y1,z1,w1)

F1 =Matrix([0,(a*c)/(b*(c+a)),c*t/a,0])

F = Matrix([  d * y -b*x*(z+w),
    a*(x+z)- b*y*(z+w) -(c+d)*x,
    
    b*x*(z+w) + d*w -e*z,
  b*y*(z+w) - (c+e+d)*w] )

 
A = ff.Joc(F, Var, 4)

B = A.subs({x:F1[0], y:F1[1], z: F1[2], w: F1[3]})

B1 = B.subs({d:0,e:0})

ev =B1.eigenvals()

evv =B1.eigenvects()

q1= Matrix([0,0,c/a, 1])

BT = B1.transpose()

evv2 = BT.eigenvects()

p1 = Matrix([1,0,1,0])
 

p= (a/c)*p1

pq =p.dot(q1)




MB=ff.MulB(F,Var,Var,4)

#Binv=B.inv()
 
MB1 = (MB.subs({x:q1[0],y:q1[1],z:q1[2], w: q1[3]}))


