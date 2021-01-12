clear;clc;close all;

%% options
makegif=0;
showanimation=1;
kymograph=0;
drawperframe=100;
L=100; % half-domain size
nx=200;
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
Wmax = 1.0; %1.0
%W=@(x,t) Wmax*heaviside(x-(-80+growthrate*max(t-50,0)));
W=@(x,t) ones(size(x))*0;

uyy_est = -0.0; % -0.28810;
vyy_est =  0.0; %  0.06471;

%% reaction
f = @(u,v,x,t) gamma * (a + W(x,t) - u + (u.^2).*v) + Du*uyy_est;
g = @(u,v,x,t) gamma * (b - (u.^2).*v) + Dv*vyy_est;
u0 = a+Wmax+b+Du\gamma*uyy_est+Dv\gamma*vyy_est;
v0 = (b+Dv\gamma*vyy_est)/(u0^2);
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

if ispc % is windows
    folder='D:\liuyueFolderOxford1\turingpattern\simulations\';
else % is linux
    folder='/home/liuy1/Documents/turingpattern/simulations/';
end
if uyy_est==0 && vyy_est==0
    uyytext = '';
else
    uyytext = 'uyyAdjusted_';
end
prefix = strcat('schnackenberg_1d_',uyytext , datestr(datetime('now'), 'yyyymmdd_HHMMSS'),'_a=',num2str(a), '_b=', num2str(b), '_growth=', num2str(growthrate) );
prefix = strcat(folder, prefix);
if makegif
    save([prefix,'.mat'], '-mat');
end

%% Set up figure
giffile = [prefix,'.gif'];
if showanimation
    fig_pos = [10 100 1000 500];
    fig=figure('Position',fig_pos,'color','w');
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
    tightEdge(gca);
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

%% estimate period and amplitude
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

[peakval,peakloc]=findpeaks(u); %does not include possible peaks on boundary
uxx=u(peakloc+1)+u(peakloc-1)-2*u(peakloc);
vxx=v(peakloc+1)+v(peakloc-1)-2*v(peakloc);
uxx_peaks = mean(uxx);
vxx_peaks = mean(vxx);
fprintf('u_{xx} Averaged on peaks: %.5f\n',uxx_peaks);
fprintf('v_{xx} Averaged on peaks: %.5f\n',vxx_peaks);

[valyval,valyloc]=findpeaks(-u); %does not include possible peaks on boundary
uxx=u(valyloc+1)+u(valyloc-1)-2*u(valyloc);
vxx=v(valyloc+1)+v(valyloc-1)-2*v(valyloc);
uxx_valy = mean(uxx);
vxx_valy = mean(vxx);
fprintf('u_{xx} Averaged on valleys: %.5f\n',uxx_valy);
fprintf('v_{xx} Averaged on valleys: %.5f\n',vxx_valy);

%% save
if makegif
    ufinal = u;
    vfinal = v;
    save([prefix,'.mat'],'ufinal','vfinal','uavg','uamp','vavg','vamp','est_period','uxx_peaks','vxx_peaks','uxx_valy','vxx_valy', '-mat','-append');
end
if kymograph
    kymograph_pos = [100,100,650,500];
    u_kymograph = plot_kymograph(uu, kymograph_pos,T,[-L,L],'u');
    v_kymograph = plot_kymograph(vv, kymograph_pos,T,[-L,L],'v');
    saveas(u_kymograph,[prefix,'_ukymograph.png']);
    saveas(v_kymograph,[prefix,'_vkymograph.png']);
end
