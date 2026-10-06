from fractions import Fraction as F
from math import gcd
from functools import reduce
# qBL = 3(B-L), r = 2 T3R, y = 6Y = qBL + 3r
states=[
 ("Q_up",1,0,F(1,2),1, F(2,3)),
 ("Q_down",1,0,F(-1,2),1,F(-1,3)),
 ("L_nu",-3,0,F(1,2),-3,F(0)),
 ("L_e",-3,0,F(-1,2),-3,F(-1)),
 ("u^c",-1,-1,F(0),-4,F(-2,3)),
 ("d^c",-1,1,F(0),2,F(1,3)),
 ("e^c",3,1,F(0),6,F(1)),
 ("nu^c",3,-1,F(0),0,F(0))]
ys=[];qbs=[]
for name,q,r,t3,y,qem in states:
 assert q+3*r==y
 assert t3+F(y,6)==qem
 ys.append(abs(y));qbs.append(abs(q))
 print(name,'qBL',q,'r',r,'6Y',y,'Qem',qem)
assert reduce(gcd,[x for x in ys if x])==1
assert reduce(gcd,[x for x in qbs if x])==1
print('primitive qBL gcd',reduce(gcd,[x for x in qbs if x]))
print('primitive y gcd',reduce(gcd,[x for x in ys if x]))
print('hypercharge quantum',F(1,6))
print('PASS')
