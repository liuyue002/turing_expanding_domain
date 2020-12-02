clear;clc;close all;

%% options
makegif=1;
showanimation=1;
drawperframe=500;
L=100; % half-domain size
nx=200;
dx=2*L/nx;
growthrate = 0.1; % bif, 0.05 to 0.75
T=200/growthrate + 50;
%T=100;
dt=0.002;
nt=T/dt+1;
nFrame=ceil((T/dt)/drawperframe);

%% parameters
gamma = 5;
a = 0.2;
b = 2.0;
Du = 1;
Dv = 100;
%growthrate = 0.2; % bif, 0.05 to 0.75
Wmax = 0.6; %1.5
%radius = @(t) min(growthrate*t + heaviside(t-350)*0.3*(t-350), 100);
%radius = @(t) min(growthrate*t, 100);
%radius = @(t) 50 + min((growthrate * (max(t-50,0))), 50);
%radius = @(t) floor(min(growthrate * t, 100)/2.5)*2.5;
%W=@(x,y,t) ((x.^2 + y.^2) > (radius(t)^2)) * Wmax; % illumination level
%W=@(x,y,t) Wmax - (y<80).*(y> -80).*(x>-80).*(x < (-60 + growthrate*t)) *Wmax;
W=@(x,y,t) Wmax*heaviside(x-(-80+growthrate*t));
%W=@(x,y,t) ones(size(x))*0;

%% reaction
f = @(u,v,x,y,t) gamma * (a + W(x,y,t) - u + (u.^2).*v);
g = @(u,v,x,y,t) gamma * (b - (u.^2).*v);
u0 = a+Wmax+b;
v0 = b/(u0^2);
noisestrength = 0.01;
fprintf('Equilibrium: u0=%.5f, v0=%.5f\n',u0,v0);

%% FDM setup
x=linspace(-L,L,nx)';
y=linspace(-L,L,nx)';
[X,Y] = meshgrid(x,y);
xmeshvec = reshape(X,[nx^2,1]);
ymeshvec = reshape(Y,[nx^2,1]);
u=zeros(nx,nx);
v=zeros(nx,nx);

I=speye(nx);
e=ones(nx,1);
T = spdiags([e -4*e e],[-1 0 1],nx,nx);
T(1,1)=-3;
T(end,end)=-3;
S = spdiags([e e],[-1 1],nx,nx);
A = (kron(I,T) + kron(S,I));
T2=spdiags([e,-3*e,e],[-1 0 1],nx,nx);
T2(1,1)=-2;
T2(nx,nx)=-2;
A(1:nx,1:nx)=T2;
A(end-nx+1:end,end-nx+1:end)=T2;
A = A/(dx^2);

%% initial condition
u(:)=u0;
%u = u + (rand(size(u))*0.6-0.3);
%u = rand(size(u))*3;
q=0.6;
u = 2.5 + 1.5*cos(q*Y);
v(:)=v0;

folder='/home/liuy1/Documents/turingpattern/simulations/';
prefix = strcat('schnackenburg_2d_noisy_square_wavyinit_' , datestr(datetime('now'), 'yyyymmdd_HHMMSS'), '_b=', num2str(b), '_growth=', num2str(growthrate) );
prefix = strcat(folder, prefix);
if makegif
    save([prefix,'.mat'], '-mat');
end

%% Set up figure
giffile = [prefix,'.gif'];
if showanimation
    fig_pos = [10 400 1660 550];
    fig=figure('Position',fig_pos);
    numticks=10;
    sfig1=subplot(1,3,1);
    ufig=imagesc(u,[0,6]);
    set(gca,'YDir','normal');
    colormap('hot');
    colorbar;
    axis image;
    set(sfig1,'XTick',0:(nx/numticks):nx);
    set(sfig1,'YTick',0:(nx/numticks):nx);
    set(sfig1,'XTickLabel',num2str((-L:2*L/numticks:L)'));
    set(sfig1,'YTickLabel',num2str((-L:2*L/numticks:L)'));
    xlim([0,200]);
    ylim([0,200]);
    utitle=title('u, t=0');
    pbaspect([1 1 1]);
    %ucirc=draw_circle(nx/2,nx/2,0);
    sfig2=subplot(1,3,2);
    vfig=imagesc(v, [0, 0.8]);
    set(gca,'YDir','normal');
    colormap('hot');
    colorbar;
    axis image;
    set(sfig2,'XTick',0:(nx/numticks):nx);
    set(sfig2,'YTick',0:(nx/numticks):nx);
    set(sfig2,'XTickLabel',num2str((-L:2*L/numticks:L)'));
    set(sfig2,'YTickLabel',num2str((-L:2*L/numticks:L)'));
    xlim([0,200]);
    ylim([0,200]);
    title('v');
    pbaspect([1 1 1]);
    %vcirc=draw_circle(nx/2,nx/2,0);
    sfig3=subplot(1,3,3);
    Wval=W(X,Y,0);
    Wfig=imagesc(Wval,[0, Wmax]);
    set(gca,'YDir','normal');
    colormap('hot');
    colorbar;
    axis image;
    set(sfig3,'XTick',0:(nx/numticks):nx);
    set(sfig3,'YTick',0:(nx/numticks):nx);
    set(sfig3,'XTickLabel',num2str((-L:2*L/numticks:L)'));
    set(sfig3,'YTickLabel',num2str((-L:2*L/numticks:L)'));
    xlim([0,200]);
    ylim([0,200]);
    title('W');
    pbaspect([1 1 1]);
end

%% simulation iteration
th=0.5; % 0: fw euler, 0.5: Crank-Nicosen, 1: bw euler
Tu=speye(nx^2)-th*dt*Du*A;
Tv=speye(nx^2)-th*dt*Dv*A;

for ti=1:1:nt
    t=dt*(ti-1);
    if (mod(ti, drawperframe) == 1)
        if showanimation
            ufig.CData=u;
            vfig.CData=v;
            utitle.String=['u, t=',num2str(t)];
            %update_circle(ucirc,nx/2,nx/2,radius(t) * (nx/(2*L)));
            %update_circle(vcirc,nx/2,nx/2,radius(t) * (nx/(2*L)));
            Wval=W(X,Y,t);
            Wfig.CData=Wval;
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
    end
    
    uvec=reshape(u,[nx^2,1]);
    vvec=reshape(v,[nx^2,1]);
    fvec=f(uvec,vvec,xmeshvec,ymeshvec,t);
    gvec=g(uvec,vvec,xmeshvec,ymeshvec,t);
    
    urhs = uvec + dt*(fvec + (1-th)*Du*A*uvec);
    unew = Tu\urhs;
    u = reshape(unew,[nx,nx]);
    u = u + normrnd(0,noisestrength,size(u));
    
    vrhs = vvec + dt*(gvec + (1-th)*Dv*A*vvec);
    vnew = Tv\vrhs;
    v = reshape(vnew,[nx,nx]);
    
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

