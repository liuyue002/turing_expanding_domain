% calculate u_{yy} at the croos section at y=0 (horizontal cross section)

%% estimate from 2D simulations, with noise in u 
load('simulations/schnackenburg_2d_noisy_narrow_hssinit_20201130_144415_a=0.05_b=1.4_growth=0.1.mat');
u=ufinal;
v=vfinal;

uyy = u(round(nx/2)+1,:) + u(round(nx/2)-1,:) -2*u(round(nx/2),:);
vyy = v(round(nx/2)+1,:) + v(round(nx/2)-1,:) -2*v(round(nx/2),:);
figure;
hold on;
plot(x,uyy);
plot(x,vyy);
legend('u_{yy}', 'v_{yy}');
xlabel('x');
hold off;
fprintf('u_{yy} Average: %.5f\n',mean(uyy));
fprintf('v_{yy} Average: %.5f\n',mean(vyy));

%% estimate from 1D sim
load('simulations/schnackenburg_1d_20201202_130108_b=1.4_growth=0.1.mat')
u=ufinal;
v=vfinal;

fig_pos = [10 100 1300 500];
figure('Position',fig_pos);
hold on
xlabel('x');
ylabel('u,v');
axis([-L,L,0,4]);
ufig=plot(x,u);
vfig=plot(x,v);
legend('u','v');
hold off

[peakval,peakloc]=findpeaks(u); %does not include possible peaks on boundary
uyy=u(peakloc+1)+u(peakloc-1)-2*u(peakloc);
vyy=v(peakloc+1)+v(peakloc-1)-2*v(peakloc);
fprintf('u_{yy} Average: %.5f\n',mean(uyy));
fprintf('v_{yy} Average: %.5f\n',mean(vyy));