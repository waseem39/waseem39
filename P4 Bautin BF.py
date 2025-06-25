# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 16:35:03 2024

@author: Wasim
"""
from sympy import*

import cmfff as ff

x, y, z, w, x1, y1, z1, w1, = symbols('x y z w x1 y1 z1 w1')

a, b, c , d , e ,k, W = symbols('a, b, c, d, e ,k W')
 
Var = (x,y,z,w)

Var1 = (x1,y1,z1,w1)


F1 =Matrix([k,0,0,0])

F = Matrix([  d * y -b*x*(z+w),
    a*(x+z)- b*y*(z+w) -(c+d)*x,
    b*x*(z+w) + d*w -e*z,
  b*y*(z+w) - (c+e+d)*w] )

 
A = ff.Joc(F, Var, 4)

B1 = A.subs({x:F1[0], y:F1[1], z: F1[2], w: F1[3]})


B2 = B1.subs({a:0,b:0})

B = B2.applyfunc(simplify)


ev = B.eigenvals()

evv =B.eigenvects()


q1 = Matrix([W*I, 1,0, 0]) # sqrt(-d*(c + d))/(c + d)=Wi


q2 = Matrix([-W*I, 1,0, 0])

BT = B.transpose()

ev1 = BT.eigenvals()

evv2 =BT.eigenvects()

p = Matrix([-W*I, 1,0, 0])

p1 = (1/(1-W**2 ))*p

p2 = Matrix([W*I, 1,0, 0])

pq =p1.dot(q1) #to veufy simplify(pq) 

MB=ff.MulB(F,Var,Var1,4)

Binv=B.inv()
 
MB1 = (MB.subs({x:q1[0],y:q1[1],z:q1[2], w: q1[3],x1:q2[0],y1:q2[1],z1:q2[2], w1:q2[3] }))

s = Binv* MB1 

s1= MB.subs({x:q1[0], y:q1[1],z:q1[2], w:q1[3], x1:s[0], y1:s[1], z1: s[2],w1:s[3]})

s2 = s1.subs({a :0,b:0})

l1 =p2.dot(s2)


r = ((2*I*W*eye(4) - B).inv())


MB2 = (MB.subs({x:q1[0],y:q1[1],z:q1[2], w : q1[3] , x1:q1[0],y1:q1[1],z1:q1[2],w1: q1[3]}))

r1 = r*MB2

r3 = MB.subs({x:q2[0], y:q2[1],z:q2[2], w:q2[3],x1:r1[0], y1:r1[1], z1: r1[2],w1:r1[3] })


l2 = p1.dot(r3)


l3 = simplify(l2)

lye = re(l3 -2*l1)

lyex =    simplify(lye)


lyexx = (1/(2*W))*(lyex)


LEC = simplify(lyexx)   

# For 2nd lyapunove coefficent we proceed as

# we have trilinear and higher multiliner function are zero

#And B(u,u) is zero for all values hence 2nd coefficent are zero we have generic case