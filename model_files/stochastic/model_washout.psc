# Reactions
R1:
    p_off > p_on
    k_on*p_off

R2:
    p_on > p_off
    k_off*p_on

R3:
    p_on > c + p_on
    ksyn_c20*p_on
R41:
    p_off2 > p_on2
    k_on2*p_off2

R42:
    p_on2 > p_off2
    k_off2*p_on2

R43:
    p_on2 > x + p_on2
    ksyn_x*p_on2

	
R5:
    c > $pool
    kdegbg*c
	
R6:
	mstar > m
	kinact*mstar/(J+mstar)
	
	
R7:
	m  > mstar
	kact*x*m*k_nuk*nuk/((J1+nuk)*(J+m))

R8:
    mstar + c > mc
    kassmc*mstar*c
	
R9:
    mc > mstar + c
    kdissmc*mc

R10:
    a + c > ac
    kassac*a*c

R11:
    mc + ac > acmc
    kassacmc*ac*mc

R12:
    ac > a + c
    kdissac*ac

R13:
    acmc > mc + ac
    kdissacmc*acmc

R14:
    acmc > a + c + m
    kdeg*acmc

R15:
    x + ac > ac
    kdeg_x*x*ac
	
R16:
	mc > m
	kdegbg*mc
	
R17:
	ac > a
	kdegbg*ac 

R18:
	acmc > a + m
	kdegbg*acmc
	
R19:
	x > $pool
	kdegbg_x*x
R20:
	nuk > $pool 
	nuk*knUK_att*P
	

# InitPar

ksyn_c20 = 46.9619
ksyn_x = 12.5
kdegbg = 0.03808
kinact = 0.0067
J = 2.1127
kact = 0.0242
k_nuk = 169.9297
J1 = 0.4
kassmc = 0.0215
kdissmc = 0.0313
kassac = 0.1073
kassacmc = 0.5794
kdissac = 0.9414
kdissacmc = 0.6115
kdeg = 1.3127
kdeg_x = 0.024
kdegbg_x = 0.0617

nuk=10



k_on = 1
k_off = 0.1
k_on2 = 1
k_off2 = 0.1

knUK_att=0.05
P=0
# InitVar
p_off = 0
p_on = 1
p_off2 = 0
p_on2 = 1

c = 13
mstar = 140
mc = 4
ac = 25
acmc = 29 
m = 2
a = 45
x = 17


# Assignment rules
!F ctot = c + mc + ac + 2*acmc
P=0
FIX: P


Event: start_attachment, _TIME_ >=20.0, 0.0


{
P=1
}

