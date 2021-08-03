%% options
% file path where the results are saved
prefix='~/growing_domain_simulation';
drawperframe=50;
L=100; % half-domain size
nx=400;
dx=2*L/nx;
growthrate = 0; % bif, 0.05 to 0.75
if growthrate == 0
    T=200;
else
    T=200/growthrate + 100;
end
dt=0.01;
nt=T/dt+1;
nFrame=ceil((T/dt)/drawperframe);
sympref('HeavisideAtOrigin',0);

%% parameters
a=12;
b=0.31;
d=1;
sigma=50;
Wmax = 1.5;
if growthrate == 0
    rho=@(t) L+10; % put L+10 for full-sized domain to avoid annoying boundary issue
    W=@(x,t) ones(size(x))*0;
    Wmax=0;
else
    rho=@(t) -L+growthrate*t;
    W=@(x,t) Wmax*heaviside(x-rho(t));
end

%% reaction
f = @(u,v,x,t) a-u-4*u.*v./(1+u.^2)-W(x,t);
g = @(u,v,x,t) sigma*b*(u-u.*v./(1+u.^2) + W(x,t));
u0 = (a-5*Wmax)/5;
v0 = (u0+Wmax)*(1+u0^2)/u0;
noisestrength = 0; % 0 or 0.01
fprintf('Equilibrium outside of effective domain: u0=%.5f, v0=%.5f\n',u0,v0);

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
u(1)=2*u0;
v(:)=v0;

save([prefix,'.mat'], '-mat');

%% Set up figure
giffile = [prefix,'.gif'];
fig_pos = [10 100 1000 500];
fig=figure('Position',fig_pos,'color','w');
hold on
xlabel('x');
ylabel('u,v,W');
axis([-L,L,0,8]);
ufig=plot(x,u);
vfig=plot(x,v);
Wval = W(x,0);
wfig=plot(x,Wval);
legend('u','v','W');
figtitle=title('t=0');
tightEdge(gca);
hold off

uu=zeros(nFrame,nx);
vv=zeros(nFrame,nx);

%% simulation iteration
th=0.5; % 0: fw euler, 0.5: Crank-Nicosen, 1: bw euler
Du = 1;
Dv = d*sigma;
Tu=speye(nx)-th*dt*Du*A;
Tv=speye(nx)-th*dt*Dv*A;

for ti=1:1:nt
    t=dt*(ti-1);
    if (mod(ti, drawperframe) == 1)
        ufig.YData=u;
        vfig.YData=v;
        Wval=W(x,t);
        wfig.YData=Wval;
        figtitle.String=['t=',num2str(t,'%.1f')];
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
        
        uu(iFrame,:)=u;
        vv(iFrame,:)=v;
    end
    
    fvec=f(u,v,x,t);
    gvec=g(u,v,x,t);
    
    urhs = u + dt*(fvec + (1-th)*Du*A*u);
    unew = Tu\urhs;
    unew = unew + normrnd(0,noisestrength,size(u));
    
    vrhs = v + dt*(gvec + (1-th)*Dv*A*v);
    vnew = Tv\vrhs;
    
    u=unew; v=vnew;
end

%% plot kymographs and save
ufinal = u;
vfinal = v;
save([prefix,'.mat'],'ufinal','vfinal','-mat','-append');

kymograph_pos = [100,100,650,500];
u_kymograph = plot_kymograph(uu, kymograph_pos,T,[-L,L],NaN,'u',0);
v_kymograph = plot_kymograph(vv, kymograph_pos,T,[-L,L],NaN,'v',0);
saveas(u_kymograph,[prefix,'_ukymograph.png']);
saveas(v_kymograph,[prefix,'_vkymograph.png']);
