R1:
    APCi > APCa
    kap*APCi
    
R2:
    APCi > APCa
    ka*(APCa^n)/(J^n+APCa^n)*APCi
    
R3:
    APCa > APCi
    ki*nuk*APCa
    

 
#parameters    
J=27.5
ka=1
ki=0.3
kap=0.005
nuk=4
APCtot=100
n=3

#init cond
# cannot find the steady state
CycB=20
APCa=1
APCi=99

#assignment rules
APCi = APCtot - APCa

P=0
FIX: P
Event: start_attachment, _TIME_ >=100.0, 0.0
{
P=1
}
