#!/usr/bin/env python3

N=20+1
c=['']*N

for n in range(0,N):
    cn=''
    for i in range(0,n+1):
        for j in range(0,n-i+1):
                cn+='u'+str(i)+'*u'+str(j)+'*v'+str(n-i-j)+'+'
    cn=cn[:-1]
    c[n]=cn

with open('schnackenberg_fourier.txt','w') as f:
    print('c0='+c[0],file=f)
    print('u0\'=a-('+c[0]+')',file=f)
    print('v0\'=b-v0+'+c[0],file=f)
    for n in range(1,N):
        print('c'+str(n)+'='+c[n],file=f)
        ueqn='u'+str(n)+'\'='+'-u'+str(n)+'*'+str(n**2)+'*(pi^2)*du'+'/(L^2)'+'-c'+str(n)
        print(ueqn,file=f)
        veqn='v'+str(n)+'\'='+'-v'+str(n)+'*('+str(n**2)+'*(pi^2)*dv'+'/(L^2) +1)+'+'c'+str(n)
        print(veqn,file=f)
#    vbd='vbd='
#    for n in range(0,N):
#        vbd += 'v'+str(n)+'+'
#    vbd=vbd[:-1]
#    print(vbd,file=f)
        
with open('schnackenberg_fourier.txt','r') as f:
    kinetics=f.read()
        
with open('schnackenberg_template.txt','r') as f2:
    template=f2.read()
    
with open('schnackenberg_fourier.ode','w') as f3:
    print(template.replace('${KINETICS}',kinetics).replace('${N}',str(N-1)),file=f3)