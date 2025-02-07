R1:
    APCi > APCa
    ka*APCi
    
R2:
    APCa > APCi
    ki*X*nuk*APCa
    
R3:
    $pool > X
    ks
    
R4: 
    X > $pool
    kdeg*APCa^n/(J^n+APCa^n)*X
    
R5:
    X > $pool
    kdegp*X
 
#parameters    
J=10
ka=0.3
ki=0.2
nuk=4
APCtot=100
ks=0.08
kdeg=0.6
kdegp=0.02
n=3

#init cond
X=20
APCa=1
APCi=99

#assignment rules
APCi = APCtot - APCa

P=0
FIX: P

Event: start_attachment, _TIME_ >=5000.0, 0.0


{
P=1
}

