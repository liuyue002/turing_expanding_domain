%% well-mixed system analysis
syms u v q q2

% my params
gamma = 1;a = 0.05;b = 1.4;du = 1;dv = 20;
% Barrass params
%gamma = 1;a = 0.1;b = 0.9;du = 0.06;dv = 1;

W=0.0;
alpha=0;
uyy_est = -0.28810;
vyy_est =  0.06471;
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
qval = linspace(0,2,500);
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

%% find where eigenvalue crosses 0
[~,ind]=findpeaks(-abs(eigenreal));
if size(ind,2)==2
    qbif1=qval(ind(1));
    qbif2=qval(ind(2));
    fprintf('Unstable between q=%.3f and q=%.3f\n',qbif1,qbif2);
    fprintf('If quantized: q=n*pi/L\n');
    L=1;
    fprintf('If L=%.2f, this corresponds to n=%.3f to n=%.3f\n',L,qbif1*L/pi,qbif2*L/pi);
end

%% travelling wave stuff
syms c
Jtw = [0,0,1,0;
       0,0,0,1;
       1-2*u*v, -u^2, -c/du,0;
       2*u*v, u^2, 0, -c/dv];