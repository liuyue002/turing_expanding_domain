%% well-mixed system analysis
syms u v q q2
F=0.0295;
K=0.0561;
du = 0.16;
dv = 0.08;
W=0.0;
f = -u*v^2 + F*(1-u)-W;
g =  u*v^2 -(F+K)*v;
d=1-(4*(F+K)^2)/F;
u0 = 1;
v0 = 0;

Jwm = [diff(f,u), diff(f,v);
       diff(g,u), diff(g,v)];
   
disp(eig(double(subs(Jwm, [u,v], [u0,v0]))));

%% Turing analysis
D = [du,0; 0, dv];
M = Jwm-D*q2;
detM = det(subs(M, [u,v], [u0,v0]));
coefs = coeffs(diff(detM, q2));
q2min = -coefs(1)/coefs(2); % q2min where the det is at minimum
disp(double(q2min));

% q2val = linspace(-2,20,100);
% fig=figure;
% plot(q2val, subs(detM, q2, q2val));
% hline=refline(0,0);
% hline.Color = 'k';
% xlabel('q^2');
% ylabel('detM');
% ylim([-30,30]); % if det>0, then turing stable

%% plot real part of dominant eigenvalue against q
M2 = Jwm-D*q^2;
qval = linspace(0,5,50);
eigM = eig(subs(M2, [u,v], [u0,v0]));
eigenreal = max(real(double(subs(eigM,q,qval))),[],1);
fig=figure;
plot(qval, eigenreal);
rline=refline(0,0);
rline.Color=[0,0,0];
xlabel('q (wave number)');
ylabel('Re(\lambda)');
ylim([-1,1]);

[eigenmax,maxind] = max(eigenreal);
fprintf('Max eigenvalue occurs at q = %.3f, corresponding period %.3f\n', qval(maxind), 2*pi/qval(maxind));