clear;clc;close all;

%% options
makegif=0;
showanimation=1;
kymograph=0;
drawperframe=500;
L=100; % half-domain size
nx=200;
dx=2*L/nx;
%growthrate = 0.1; % bif, 0.05 to 0.75
%T=200/growthrate + 50;
T=200;
dt=0.001;
nt=T/dt+1;
nFrame=ceil((T/dt)/drawperframe);

%% parameters
F=0.0392;
K=0.0649;
Du = 0.16;
Dv = 0.08;
growthrate = 0.0; % bif, 0.05 to 0.75
Wmax = 1.0; %1.5
%radius = @(t) min(growthrate*t + heaviside(t-350)*0.3*(t-350), 100);
%radius = @(t) min(growthrate*t, 100);
%radius = @(t) 50 + min((growthrate * (max(t-50,0))), 50);
%radius = @(t) floor(min(growthrate * t, 100)/2.5)*2.5;
%W=@(x,y,t) ((x.^2 + y.^2) > (radius(t)^2)) * Wmax; % illumination level
%W=@(x,y,t) Wmax - (y<80).*(y> -80).*(x>-80).*(x < (-60 + growthrate*t)) *Wmax;
%W=@(x,y,t) Wmax*heaviside(x-(-80+growthrate*t));
W=@(x,y,t) ones(size(x))*0;

%% reaction
f = @(u,v,x,t) -(v.^2).*u + F*(1-u) - W(x,t);
g = @(u,v,x,t)  (v.^2).*u - (F+K)*v;
d=1-(4*(F+K)^2)/F;
if d>= 0 
    u0 = (1-sqrt(d))/2;
    v0 = (F/(F+K)*(1+sqrt(d)))/2;
else
    u0=1;
    v0=1;
    fprintf('Warning: nontrivial HSS do not exist\n');
end
u0=1;v0=0;
noisestrength = 0.0;
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
%u(:)=u0;
u = u0*(rand(size(u))*0.6+0.7);
%u = rand(size(u))*3;
%q=2*pi*0.09;
%u = 1.64 + 0.70*cos(q*Y);
%v(:)=v0;
%v = v0*(rand(size(v))*0.6+0.7);
v = rand(size(v))*0.4;
%v = 0.60 - 0.13*cos(q*Y);

if ispc % is windows
    folder='D:\liuyueFolderOxford1\turingpattern\simulations\';
else % is linux
    folder='/home/liuy1/Documents/turingpattern/simulations/';
end
prefix = strcat('grayscott_1d_noisy_square_hssinit_' , datestr(datetime('now'), 'yyyymmdd_HHMMSS'),'_F=',num2str(F), '_K=', num2str(K), '_growth=', num2str(growthrate) );
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
    axis([-L,L,0,3]);
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
            utitle.String=['u, t=',num2str(t)];
            Wval=W(x,t);
            Wfig.YData=Wval;
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
    u = u + normrnd(0,noisestrength,size(u));
    
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
