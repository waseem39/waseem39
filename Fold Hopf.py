# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 12:24:05 2024

@author: Wasim
"""

from sympy import*

import cmfff as ff

x, y, z, w, x1, y1, z1, w1 = symbols('x y z w x1 y1 z1 w1')

a, b, c , d , e ,k, u, v, W   = symbols('a, b, c, d, e ,k, u, v,  W')
 
Var = (x,y,z,w)

Var1 = (x1,y1,z1,w1)
F1 =Matrix([0,0,0,0])

F = Matrix([  d * y -b*x*(z+w),
    a*(x+z)- b*y*(z+w) -(c+d)*x,
    b*x*(z+w) + d*w -e*z,
  b*y*(z+w) - (c+e+d)*w] )

 
A = ff.Joc(F, Var, 4)

B = A.subs({x:F1[0], y:F1[1], z: F1[2], w: F1[3]})

B1 = B.subs({b:0,e:0})

ev =B1.eigenvals()

evv =B1.eigenvects()


#Adjoint eigen vector 


BT = B1.transpose()

ev2 = BT.eigenvals()

evv2=BT.eigenvects()



q0 = Matrix([a / W ,0,1,0])

q1 = Matrix([-(sqrt(d)*I)/ W, 1,0, 0]) # sqrt(-d*(c + d))/(c + d)=Wi

q2 =  Matrix([(sqrt(d)*I)/ W, 1,0, 0])

# #Adjoint eigeb vectors


p0 = Matrix([0,0,1, d /(c + d)])

p = Matrix([(sqrt(d)*I*W**(3) + (c+d)**2 - a*(c+d)) / (a*d),-W**2 /a -(sqrt(d)*I*W*(c+d))/ (a*d)  ,(sqrt(d)*I*W + (c+d))/d , 1])

D = (W*a*sqrt(d)) / (I*(-W**2*c - W**2*d + a*c + a*d - c**2 - 2*c*d - d**2))

P = Matrix([(-sqrt(d)*I*W**(3) - (c+d)**2 +a*(c+d)) / (a*d),-W**2 /a -(sqrt(d)*I*W*(c+d))/ (a*d)  ,(sqrt(d)*I*W + (c+d))/d , 1])


p1 = D*p

p2 = D*P

pq = p1.dot(q1)

#simplify(pq)
 # impotant to note we have to find dodt product in Complex so we
# #we have to change sign of complex values
#Bin = B1.inv()

MB=ff.MulB(F,Var,Var1,4)

#B(q0,q0)
MB1=MB.subs({x:q0[0],y:q0[1],z:q0[2],w:q0[3], x1:q0[0],y1:q0[1],z1:q0[2],w1:q0[3]})

#B(q0,q1)
MB2=MB.subs({x:q0[0],y:q0[1],z:q0[2],w:q0[3], x1:q1[0],y1:q1[1],z1:q1[2],w1:q1[3]})

#B(q1,q2)
MB3=MB.subs({x:q1[0],y:q1[1],z:q1[2],w:q1[3], x1:q2[0],y1:q2[1],z1:q2[2],w1:q2[3]})


# in our case h200 and h011 and h020 are zero
#to compute h110

C200 =( p0.dot(MB1))/2

C011 =  p0.dot(MB3)

H110 = p1.dot(MB2)

# rr = (I*W*eye(4)-B1).inv()


# r = rr.applyfunc(simplify)


# rr1 = (MB2 - p1.dot(MB1)*q1)


# r1 = rr1.applyfunc(simplify)


# h110 = (r*r1)


# hhhh = h110.applyfunc(simplify)
# h110b = h110.conjugate

# h200 = Matrix([0,0,0,0])

# h020 = Matrix([0,0,0,0])


# h011 = Matrix([0,0,0,0])

# # aaa=((p.dot(MB2))*q0-MB1)

# # Var2 = Matrix([x,I*W*x,z,w])

# # B2=B1*Var2
# # # B2 = B3.applyfunc(simplify)

# # ss=solve((Eq(B2[0],aaa[0]),Eq(B2[1],aaa[1]),Eq(B2[2],aaa[2]),Eq(B2[3],aaa[3])),(x,y,z,w))




# G200  = (p.dot(MB1))/2

# H110 = simplify(p1.dot(MB2))

# G300 = 0


# MB4=2*MB.subs({x:q0[0],y:q0[1],z:q0[2],w:q0[3], x1:h110[0],y1:h110[1],z1:h110[2],w1:h110[3]})

# H210 =( p1.dot(MB4))/2




# MB5=2*MB.subs({x:q1[0],y:q1[1],z:q1[2],w:q1[3], x1:h110[0],y1:h110[1],z1:h110[2],w1:h110[3]})


# H021 =q1.dot(MB5)

# MB6 =MB.subs({x:q1[0],y:q1[1],z:q1[2],w:q1[3], x1:h110b[0],y1:h110b[1],z1:h110b[2],w1:h110b[3]})
