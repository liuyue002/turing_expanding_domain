#!/usr/bin/env python3

N=20+1
c=['']*N #coefficients of u^2*v
u2coef=['']*N # coefficients of u^2

def u(n):
    return 'u'+str(n)

def v(n):
    return 'v'+str(n)

def u2(n):
    return 'usq'+str(n)
    #return '('+u2coef[n]+')'

u2coef0='u0^2+('
for i in range(1,N):
    u2coef0+=u(i)+'^2+'
u2coef0=u2coef0[:-1]+')/2'
u2coef[0]=u2coef0
    
for n in range(1,N):
    u2coefn=''
    # positive terms (correspond to cos(n+m) term)
    for i in range(0,(n-1)//2 +1):
        u2coefn+=u(i)+'*'+u(n-i)+'+'
    u2coefn=u2coefn[:-1]
    if n % 2 == 0:
        if n > 0:
            u2coefn+='+'
        u2coefn+=u(n//2)+'^2/2'
    u2coefn+='+'
    # negative terms (correspond to cos(n-m) term)
    for i in range(0,N-n):
        u2coefn+=u(i)+'*'+u(n+i)+'+'
    u2coefn='('+u2coefn[:-1]+')'
    u2coef[n]=u2coefn

c0=u2(0)+'*v0+('
for i in range(1,N):
    c0+=u2(i)+'*'+v(i)+'+'
c0=c0[:-1]+')/2'
c[0]=c0

for n in range(1,N):
    cn=''
    # positive terms (correspond to cos(n+m) term)
    for i in range(0,n+1):
        cn+=v(i)+'*'+u2(n-i)+'+'
    # negative terms (correspond to cos(n-m) term)
    for i in range(0,N-n):
        cn+=v(i)+'*'+u2(n+i)+'+'+v(n+i)+'*'+u2(i)+'+'
    cn='('+cn[:-1]+')/2'
    c[n]=cn
    

kinetics=''
for n in range(0,N):
    kinetics +='usq'+str(n)+'=('+u2coef[n]+')\n'
for n in range(0,N):
    kinetics +='c'+str(n)+'='+c[n]+'\n'
ubd='aux ubd='
for n in range(0,N):
    ubd+='u'+str(n)+'+'
ubd=ubd[:-1]
kinetics += ubd
        
with open('schnackenberg_template.txt','r') as f2:
    template=f2.read()
    
with open('schnackenberg_fourier.ode','w') as f3:
    print(template.replace('${FOURIER_TERMS}',kinetics).replace('${N}',str(N-1)),file=f3)