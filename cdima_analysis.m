%% well-mixed system analysis
syms u v b q q2
a = 12;
d = 1;
sigma = 50;
W = 1.5;
f = a-u-4*u*v/(1+u^2) -W;
g = sigma*b*(u - u*v/(1+u^2) +W);
u0 = (a-5*W)/5;
v0 = (u0+W)*(1+u0^2)/u0;

Jwm = [diff(f,u), diff(f,v);
       diff(g,u), diff(g,v)];
   
disp(eig(double(subs(Jwm, [u,v,b], [u0,v0,0.3]))));

%% Turing analysis
D = [1,0; 0, d*sigma];
M = Jwm-D*q2;
detM = det(subs(M, [u,v,b], [u0,v0,0.3]));
coefs = coeffs(diff(detM, q2));
q2min = -coefs(1)/coefs(2); % q2min where the det is at minimum
disp(q2min);

q2val = linspace(-2,2,100);
fig=figure;
plot(q2val, subs(detM, q2, q2val));
hline=refline(0,0);
hline.Color = 'k';
xlabel('q^2');
ylabel('detM');
ylim([-30,30]); % if det>0, then turing stable
title(['W=',num2str(W)]);

%% plot real part of dominant eigenvalue against q
M2 = Jwm-D*q^2;
qval = linspace(0,2,200);
eigM = eig(subs(M2, [u,v,b], [u0,v0,0.3]));
eigenreal = max(real(double(subs(eigM,q,qval))),[],1);
fig=figure;
plot(qval, eigenreal);
rline=refline(0,0);
rline.Color=[0,0,0];
xlabel('q (wave number)');
ylabel('Re(\lambda)');
ylim([-3,1]);
title(['W=',num2str(W)]);

[eigenmax,maxind] = max(eigenreal);
fprintf('Max eigenvalue occurs at q = %.3f, corresponding period %.3f\n', qval(maxind), 2*pi/qval(maxind));