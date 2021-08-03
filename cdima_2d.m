%% options
% file path where the results are saved
prefix='~/growing_domain_simulation';
drawperframe=100;
L=100; % half-domain size, the domain is [-L,L]
nx=200;
dx=2*L/nx;
growthrate = 0.1; % 0.1 to 2; set to 0 for a fixed domain simulation
if growthrate == 0
    T=100;
else
    T=200/growthrate + 100;
end
dt=0.02;
nt=T/dt+1;
sympref('HeavisideAtOrigin',0);

%% parameters
a=12;
b=0.31;
d=1;
sigma=50;
Wmax = 1.5;
wid=L+10; % vertical width of the effective domain, default L+10 for full size
if growthrate == 0
    rho=@(t) L+10;
    W=@(x,y,t) ones(size(x))*0;
    Wmax=0;
else
    rho=@(t) -L+growthrate*t;
    W=@(x,y,t) Wmax*(1-heaviside(-x+rho(t)).*heaviside(y+wid).*heaviside(-y+wid));
end

%% reaction kinetics
f = @(u,v,x,y,t) a-u-4*u.*v./(1+u.^2)-W(x,y,t);
g = @(u,v,x,y,t) sigma*b*(u-u.*v./(1+u.^2) + W(x,y,t));
u0 = (a-5*Wmax)/5;
v0 = (u0+Wmax)*(1+u0^2)/u0;
noisestrength = 0.01; % default 0.01

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
T1 = spdiags([e -4*e e],[-1 0 1],nx,nx);
T1(1,1)=-3;
T1(end,end)=-3;
S = spdiags([e e],[-1 1],nx,nx);
A = (kron(I,T1) + kron(S,I));
T2=spdiags([e,-3*e,e],[-1 0 1],nx,nx);
T2(1,1)=-2;
T2(nx,nx)=-2;
A(1:nx,1:nx)=T2;
A(end-nx+1:end,end-nx+1:end)=T2;
A = A/(dx^2);

%% initial condition
u(:)=u0;
u = u + (rand(size(u))*0.6-0.3);
v(:)=v0;

%alternatively:
% q=0.911;
% u = 2.405 + 1.156*cos(q*Y);
% v = 6.668 + 0.709*cos(q*Y);

uinit=u;
vinit=v;
save([prefix,'.mat'], '-mat');

%% Set up figure
giffile = [prefix,'.gif'];
urange=[0,4];
vrange=[4.5,8];

fig_pos = [100 100 1650 500];
fig=figure('Position',fig_pos,'color','w');
numticks=10;

sfig1=subplot('Position',[0.04,0,0.28,1]);
hold on
ufig=imagesc(u,urange);
Wline=plot([(rho(0)+L)/dx,(rho(0)+L)/dx],[0,nx],'-b');
xlabel('x');
ylabel('y');
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
hold off
pbaspect([1 1 1]);

sfig2=subplot('Position',[0.37,0,0.28,1]);
vfig=imagesc(v, vrange);
xlabel('x');
ylabel('y');
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

sfig3=subplot('Position',[0.70,0,0.28,1]);
Wval=W(X,Y,0);
Wfig=imagesc(Wval,[0, Wmax]);
xlabel('x');
ylabel('y');
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

%% simulation iteration
th=0.5; % 0: fw euler, 0.5: Crank-Nicosen, 1: bw euler
Du = 1;
Dv = d*sigma;
Tu=speye(nx^2)-th*dt*Du*A;
Tv=speye(nx^2)-th*dt*Dv*A;

for ti=1:1:nt
    t=dt*(ti-1);
    if (mod(ti, drawperframe) == 1)
        ufig.CData=u;
        vfig.CData=v;
        utitle.String=['u, t=',num2str(t)];
        Wline.XData=[(rho(t)+L)/dx,(rho(t)+L)/dx];
        Wval=W(X,Y,t);
        Wfig.CData=Wval;
        drawnow;
        
        iFrame=(ti-1)/drawperframe+1;
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if iFrame==1
            imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
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
end

%% saving final pattern

ufinal = u;
vfinal = v;
save([prefix,'.mat'],'ufinal','vfinal','-mat','-append');
