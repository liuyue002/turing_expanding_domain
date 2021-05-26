%% well-mixed system analysis
syms u v q q2
gamma = 1;
a = 0.05;
b = 1.6;
du = 1;
dv = 20;
W=0.0;
alpha=0;
%uyy_est = -0.28810;
uyy_est=0;
%vyy_est =  0.06471;
vyy_est=0;
a_eff = a + alpha*du*uyy_est/gamma +W;
b_eff = b + alpha*dv*vyy_est/gamma;
f = gamma * (a_eff -u+u^2 *v );
g = gamma * (b_eff - u^2 *v );
u0 = a_eff+b_eff;
v0 = b_eff/(u0^2);

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
qval = linspace(0,1.5,500);
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


