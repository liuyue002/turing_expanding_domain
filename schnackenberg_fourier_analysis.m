syms u v u0 v0 u1 v1 L
a=0.05;
b=1.4;
du=1;
dv=20;
ueq=a+b;
veq=b/(ueq^2);
f = a -u+u^2 *v ;
g = b - u^2 *v ;
Jwm = [diff(f,u), diff(f,v);
       diff(g,u), diff(g,v)];

A = du*diff(g,v) + dv*diff(f,u);
Aval = double(subs(A,[u,v],[ueq,veq]));
fprintf('A=%.3f\n',Aval);

Lmult1=double(subs(pi*sqrt(2*dv*du/(A+sqrt(A^2-4*dv*du*det(Jwm)))), [u,v],[ueq,veq]));
fprintf('Lmult+ =%.3f\n',Lmult1);
Lmult2=double(subs(pi*sqrt(2*dv*du/(A-sqrt(A^2-4*dv*du*det(Jwm)))), [u,v],[ueq,veq]));
fprintf('Lmult- =%.3f\n',Lmult2);

%% 

c0=u0^2*v0 + u1^2*v0/2 + u0*u1*v1;
c1=2*u0*u1*v0 + u0^2*v1 + u1^2*v1/2;
f0 = a-u0+c0;
g0 = b-c0;
f1 = -u1*du*(pi/L)^2 - u1 + c1;
g1 = -v1*dv*(pi/L)^2 - c1;
vars=[u0,v0,u1,v1];
rhs =[f0,g0,f1,g1];

J = jacobian(rhs,[vars,L]);
rhsval=double(subs(rhs,[vars,L],[1,1.0145262,0.9654520,-0.6083693,2.3721403]));
Jval=double(subs(J,[vars,L],[1,1.0145262,0.9654520,-0.6083693,2.3721403]));