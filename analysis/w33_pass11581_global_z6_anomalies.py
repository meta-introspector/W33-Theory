import cmath, math
from collections import Counter

# One Weyl-16 decomposition from the exact Spin(10) probe, y=6Y.
blocks=[
 ("Q",(3,2), 1,6,"3"),
 ("L",(1,2),-3,2,"1"),
 ("uc",(3,1),-4,3,"bar3"),
 ("dc",(3,1), 2,3,"bar3"),
 ("ec",(1,1), 6,1,"1"),
 ("nc",(1,1), 0,1,"1"),
]

# Local anomalies in integer normalization y=6Y.
grav=sum(mult*y for _,_,y,mult,_ in blocks)
u1c=sum(mult*y**3 for _,_,y,mult,_ in blocks)
su2=sum((3 if name=="Q" else 1 if name=="L" else 0)*y for name,_,y,_,_ in blocks)
su3=sum((2*y if name=="Q" else y if name in ("uc","dc") else 0) for name,_,y,_,_ in blocks)
su3cube=2-1-1
witten_doublets=3+1
print("anomalies",grav,u1c,su2,su3,su3cube,witten_doublets)
assert (grav,u1c,su2,su3,su3cube)==(0,0,0,0,0)
assert witten_doublets%2==0

# Electric charge in 6Q units: q_e = y + 3 t3 where t3=+-1 on SU2 doublets.
eq=Counter()
for name,(c,w),y,mult,rep in blocks:
 if w==2:
  base=mult//2
  eq[y+3]+=base
  eq[y-3]+=base
 else:
  eq[y]+=mult
print("6Q spectrum",dict(sorted(eq.items())))
# physical left-handed conjugate convention: {Q_u,Q_d,L_nu,L_e,uc,dc,ec,nc}
assert eq==Counter({4:3,-2:3,0:2,-6:1,-4:3,2:3,6:1})

# Global kernel of SU3 x SU2 x U1_y on these reps.
# central element = (omega^a, (-1)^b, exp(i theta)), theta=n*pi/3 due e^c.
kernel=[]
omega=cmath.exp(2j*math.pi/3)
for a in range(3):
 for b in range(2):
  for n in range(6):
   theta=n*math.pi/3
   ok=True
   for name,(c,w),y,mult,rep in blocks:
    z3=(omega**a if rep=="3" else omega**(-a) if rep=="bar3" else 1)
    z2=(-1)**b if w==2 else 1
    phase=z3*z2*cmath.exp(1j*y*theta)
    if abs(phase-1)>1e-9: ok=False; break
   if ok: kernel.append((a,b,n))
print("kernel",kernel)
assert len(kernel)==6
assert set(kernel)=={(n%3,n%2,n) for n in range(6)}
print("generator",(1,1,1),"order",6)
