clear;clc;close all;

%% options
makegif=1;
showanimation=1;
kymograph=0;
drawperframe=100;
L=100; % half-domain size
nx=400;
dx=2*L/nx;
growthrate = 0.1; % bif, 0.05 to 0.75
%T=400;
T=200/growthrate + 100;
dt=0.01;
nt=T/dt+1;
nFrame=ceil((T/dt)/drawperframe);

%% parameters
gamma = 1.0;
a = 0.05;
b = 1.4;
Du = 1;
Dv = 20;
%growthrate = 0.2; % bif, 0.05 to 0.75
Wmax = 1.0; %1.5
%W=@(x,t) Wmax*heaviside(x-(-80+growthrate*t));
W=@(x,t) Wmax*heaviside(x-(-80+growthrate*max(t-50,0)));
%W=@(x,t) ones(size(x))*0;

uyy_est = -0.28810;
vyy_est =  0.06471;

%% reaction
f = @(u,v,x,t) gamma * (a + W(x,t) - u + (u.^2).*v) + Du*uyy_est;
g = @(u,v,x,t) gamma * (b - (u.^2).*v) + Dv*vyy_est;
u0 = a+Wmax+b;
v0 = b/(u0^2);
noisestrength = 0;
fprintf('Equilibrium: u0=%.5f, v0=%.5f\n',u0,v0);

%% FDM setup
x=linspace(-L,L,nx)';
u=zeros(nx,1);
v=zeros(nx,1);

o=ones(nx,1);
A=spdiags([o -2*o o],[-1 0 1],nx,nx);
A(1,1)=-1; % for no-flux BC
A(nx,nx)=-1;
A=A/(dx^2);

%% initial condition
u(:)=u0;
u = u + (rand(size(u))*0.6-0.3);
%u = rand(size(u))*3;
%q=0.8;
%u = 1.5 + cos(q*Y);
v(:)=v0;

%folder='/home/liuy1/Documents/turingpattern/simulations/';
folder='D:\liuyueFolderOxford1\turingpattern\simulations\';
prefix = strcat('schnackenburg_1d_uyyAdjusted_' , datestr(datetime('now'), 'yyyymmdd_HHMMSS'), '_b=', num2str(b), '_growth=', num2str(growthrate) );
prefix = strcat(folder, prefix);
if makegif
    save([prefix,'.mat'], '-mat');
end

%% Set up figure
giffile = [prefix,'.gif'];
if showanimation
    fig_pos = [10 100 1300 500];
    fig=figure('Position',fig_pos);
    hold on
    xlabel('x');
    ylabel('u,v,W');
    axis([-L,L,0,4]);
    ufig=plot(x,u);
    vfig=plot(x,v);
    Wval = W(x,0);
    wfig=plot(x,Wval);
    legend('u','v','W');
    figtitle=title('t=0');
    hold off
end

if kymograph
    uu=zeros(nFrame,nx);
    vv=zeros(nFrame,nx);
end

%% simulation iteration
th=0.5; % 0: fw euler, 0.5: Crank-Nicosen, 1: bw euler
Tu=speye(nx)-th*dt*Du*A;
Tv=speye(nx)-th*dt*Dv*A;

for ti=1:1:nt
    t=dt*(ti-1);
    if (mod(ti, drawperframe) == 1)
        if showanimation
            ufig.YData=u;
            vfig.YData=v;
            Wval=W(x,t);
            wfig.YData=Wval;
            figtitle.String=['t=',num2str(t)];
            drawnow;
        end
        iFrame=(ti-1)/drawperframe+1;
        if makegif
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if iFrame==1
                imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
            end
        end
        if kymograph
            uu(iFrame,:)=u;
            vv(iFrame,:)=v;
        end
    end
    
    fvec=f(u,v,x,t);
    gvec=g(u,v,x,t);
    
    urhs = u + dt*(fvec + (1-th)*Du*A*u);
    unew = Tu\urhs;
    unew = unew + normrnd(0,noisestrength,size(u));
    
    vrhs = v + dt*(gvec + (1-th)*Dv*A*v);
    vnew = Tv\vrhs;
    
    u=unew; v=vnew;
    
    if (mod(ti, drawperframe) == 0)
        utotal=sum(sum(u))*dx;
        vtotal=sum(sum(v))*dx;
        fprintf('ti=%d done, total stuff=%.2f\n',ti,utotal+vtotal);
    end
end

if makegif
    ufinal = u;
    vfinal = v;
    save([prefix,'.mat'],'ufinal','vfinal', '-mat','-append');
end
if kymograph
    plot_kymograph(uu, fig_pos, nFrame, nx,T);
    plot_kymograph(vv, fig_pos, nFrame, nx,T);
end

%%
uavg = (max(u,[],'all')+min(u,[],'all'))/2;
uamp = max(u,[],'all') - uavg;
vavg = (max(v,[],'all')+min(v,[],'all'))/2;
vamp = max(v,[],'all') - vavg;
[peakval,peakloc]=findpeaks(u);
numpeak=size(peakloc,1);
est_period=2*L/numpeak;
fprintf('For the final pattern:\n');
fprintf('uavg=%.4f, amplitude=%.4f\n',uavg,uamp);
fprintf('vavg=%.4f, amplitude=%.4f\n',vavg,vamp);
fprintf('Estimated period=%.4f, freq=%.4f\n',est_period,1/est_period);