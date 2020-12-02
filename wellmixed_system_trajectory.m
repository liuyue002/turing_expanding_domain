a=12;
b=0.3; %bif, 0.3 to 0.34
d=1;
sigma=50;
growthrate = 0.5; % bif, 0.05 to 0.75
Wmax = 1.5;

f = @(u,v,W) a-u-4*u.*v./(1+u.^2)-W;
g = @(u,v,W) sigma*b*(u-u.*v./(1+u.^2) + W);
u0 = (a-5*Wmax)/5;
v0 = (u0+Wmax)*(1+u0^2)/u0;

odefun = @(t,X) [f(X(1),X(2),Wmax); g(X(1),X(2),Wmax)];
tspan=[0 100];
init=[u0-0.1,v0];
[t,X]=ode45(odefun,tspan,init);

fig=figure;
hold on
plot(t,X(:,1));
plot(t,X(:,2));
xlabel('t');
legend('u','v');
hold off