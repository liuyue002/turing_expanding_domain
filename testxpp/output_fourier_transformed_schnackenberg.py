#!/usr/bin/env python3

N=10+1
c=['']*N

for n in range(0,N):
    cn=''
    for i in range(0,n+1):
        cn+= 'v'+str(i)+'*('
        if n-i > 0:
            cn += '2*('
            for j in range(0,(n-i-1)//2 +1):
                cn+='u'+str(j)+'*u'+str(n-i-j)+'+'
            cn=cn[:-1]
            cn+=')'
        if (n-i) % 2 == 0:
            if n-i > 0:
                cn+='+'
            cn+='u'+str((n-i)//2)+'^2'
        cn+=')+'
    cn=cn[:-1]
    c[n]=cn

kinetics=''

kinetics +='c0='+c[0]+'\n'
#kinetics +='u0\'=a-u0+c0\n'
#kinetics +='v0\'=b-c0\n'
for n in range(1,N):
    kinetics +='c'+str(n)+'='+c[n]+'\n'
#    ueqn='u'+str(n)+'\'='+'-u'+str(n)+'*('+str(n**2)+'*(pi^2)*du'+'/(L^2) +1)+'+'c'+str(n)
#    kinetics +=ueqn+'\n'
#    veqn='v'+str(n)+'\'='+'-v'+str(n)+'*'+str(n**2)+'*(pi^2)*dv'+'/(L^2)'+'-c'+str(n)
#    kinetics +=veqn+'\n'
#ubd='aux ubd='
#for n in range(0,N):
#    ubd+='u'+str(n)+'+'
#ubd=ubd[:-1]
#kinetics += ubd
        
with open('schnackenberg_template.txt','r') as f2:
    template=f2.read()
    
with open('schnackenberg_fourier.ode','w') as f3:
    print(template.replace('${FOURIER_TERMS}',kinetics).replace('${N}',str(N-1)),file=f3)