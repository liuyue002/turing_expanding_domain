
function cdima_2dfunc(growthrate)
%% options
makegif=1;
showanimation=1;
drawperframe=200;
L=100; % half-domain size
nx=200;
dx=2*L/nx;
%growthrate = 0.1; % bif, 0.05 to 0.75
if growthrate == 0
    T=100;
else
    T=200/growthrate + 100;
end
dt=0.01;
nt=T/dt+1;
sympref('HeavisideAtOrigin',0);
rng_seed=0;
rng(rng_seed); %change random seed

%% parameters
a=12;
b=0.31; %bif, 0.3 to 0.34
d=1;
sigma=50;
%growthrate = 0.2; % bif, 0.05 to 0.75
Wmax = 1.5; %1.5
wid=L+10; % put L+10 for full-sized domain to avoid annoying boundary issue
circular=0;
if growthrate == 0
    rho=@(t) L+10;
    W=@(x,y,t) ones(size(x))*0;
    Wmax=0;
else
    if circular
        rho=@(t) min(growthrate*t,L);
        W=@(x,y,t) ((x.^2 + y.^2) > (rho(t)^2)) * Wmax;
    else
        rho=@(t) -L+growthrate*t;
        W=@(x,y,t) Wmax*(1-heaviside(-x+rho(t)).*heaviside(y+wid).*heaviside(-y+wid));
    end
end

%% reaction
f = @(u,v,x,y,t) a-u-4*u.*v./(1+u.^2)-W(x,y,t);
g = @(u,v,x,y,t) sigma*b*(u-u.*v./(1+u.^2) + W(x,y,t));
u0 = (a-5*Wmax)/5;
v0 = (u0+Wmax)*(1+u0^2)/u0;
noisestrength = 0.01;
fprintf('Equilibrium outside of effective domain: u0=%.5f, v0=%.5f\n',u0,v0);

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
u(:,1)=2*u0;
%u = u + (rand(size(u))*0.6-0.3);
%u = rand(size(u))*3;
%q=0.911;
%u = 2.405 + 1.156*cos(q*Y);
v(:)=v0;
%v = 6.668 + 0.709*cos(q*Y);

if ispc % is windows
    folder='D:\liuyueFolderOxford1\turingpattern\simulations\';
else % is linux
    folder='/home/liuy1/Documents/turingpattern/simulations/';
end
if noisestrength == 0
    noisetext='nonoise_';
else
    noisetext='noisy_';
end
if circular
    circletext='circular_';
else
    circletext='';
end
if wid >= L
    widthtext='';
else
    widthtext=['narrow_wid=', num2str(wid),'_'];
end
ictext = 'boundaryinit_'; % 'hssinit_' or 'wavyinit_' or 'hssnonoiseinit_' or 'boundaryinit_'
prefix = strcat('cdima_2d_',noisetext,widthtext,circletext,ictext, datestr(datetime('now'), 'yyyymmdd_HHMMSS'),'_growth=', num2str(growthrate) );
prefix = strcat(folder, prefix);
fprintf('saving to %s\n',prefix);
if makegif
    save([prefix,'.mat'], '-mat');
end

%% Set up figure
giffile = [prefix,'.gif'];
giffile_horiz = [prefix,'_horizcrossec.gif'];
giffile_vert = [prefix,'_vertcrossec.gif'];
urange=[0,4];
vrange=[4.5,8];

if showanimation
    fig_pos = [100 100 1650 500];
    fig=figure('Position',fig_pos,'color','w');
    numticks=10;
    
    sfig1=subplot('Position',[0.04,0,0.28,1]);
    hold on
    ufig=imagesc(u,urange);
    if circular
        ucirc=draw_circle(nx/2,nx/2,0);
    else
        Wline=plot([(rho(0)+L)/dx,(rho(0)+L)/dx],[0,nx],'-b');
    end
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
    xlim([1,nx]);
    ylim([1,nx]);
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
    xlim([1,nx]);
    ylim([1,nx]);
    title('v');
    pbaspect([1 1 1]);
    
    sfig3=subplot('Position',[0.70,0,0.28,1]);
    Wval=W(X,Y,0);
    if Wmax > 0
        Wfig=imagesc(Wval,[0, Wmax]);
    else
        Wfig=imagesc(Wval,[0,1]);
    end
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
    xlim([1,nx]);
    ylim([1,nx]);
    title('W');
    pbaspect([1 1 1]);
    
    crosssec_fig_pos=[100,100,600,500];
    fig_horiz=figure('Position',crosssec_fig_pos,'color','w');
    hold on;
    horiz_crossec=plot(x,u(round(nx/2),:));
    ylim(urange);
    fig_horiz_title=title('u(y=0,t=0)');
    xlabel('x');
    ylabel('u');
    tightEdge(gca);
    hold off
    
    fig_vert=figure('Position',crosssec_fig_pos,'color','w');
    hold on;
    vert_crossec=plot(y,u(:,round(nx/2)));
    ylim(urange);
    fig_vert_title=title('u(x=0,t=0)');
    xlabel('y');
    ylabel('u');
    tightEdge(gca);
    hold off
end

%% simulation iteration
th=0.5; % 0: fw euler, 0.5: Crank-Nicosen, 1: bw euler
Du = 1;
Dv = d*sigma;
Tu=speye(nx^2)-th*dt*Du*A;
Tv=speye(nx^2)-th*dt*Dv*A;

for ti=1:1:nt
    t=dt*(ti-1);
    if (mod(ti, drawperframe) == 1)
        if showanimation
            ufig.CData=u;
            vfig.CData=v;
            utitle.String=['u, t=',num2str(t)];
            if circular
                update_circle(ucirc,nx/2,nx/2,rho(t) * (nx/(2*L)));
            else
                Wline.XData=[(rho(t)+L)/dx,(rho(t)+L)/dx];
            end
            Wval=W(X,Y,t);
            Wfig.CData=Wval;
            
            horiz_crossec.YData=u(round(nx/2),:);
            fig_horiz_title.String=['u(y=0,t=',num2str(t),')'];
            vert_crossec.YData=u(:,round(nx/2));
            fig_vert_title.String=['u(x=0,t=',num2str(t),')'];
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
            
            frame = getframe(fig_horiz);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if iFrame==1
                imwrite(imind,cm,giffile_horiz,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,giffile_horiz,'gif','WriteMode','append','DelayTime',0);
            end
            
            frame = getframe(fig_vert);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if iFrame==1
                imwrite(imind,cm,giffile_vert,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,giffile_vert,'gif','WriteMode','append','DelayTime',0);
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

%% estimate u_yy on the cross section y=0
uyy = u(round(nx/2)+1,:) + u(round(nx/2)-1,:) -2*u(round(nx/2),:);
vyy = v(round(nx/2)+1,:) + v(round(nx/2)-1,:) -2*v(round(nx/2),:);
uyy_est = mean(uyy);
vyy_est = mean(vyy);
fprintf('u_{yy}(y=0) Average: %.5f\n',uyy_est);
fprintf('v_{yy}(y=0) Average: %.5f\n',vyy_est);


%% plotting the final pattern
if makegif
    fig_pos_square = [100 100 500 500];
    fig_pos = [100 100 610 500];
    numticks=4;
    
    ufigfinal=figure('Position',fig_pos_square,'color','w');
    imagesc(u,urange);
    xlim([1,200]);
    ylim([1,200]);
    set(gca,'YDir','normal');
    axis image;
    colormap('hot');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(ufigfinal,[prefix,'_ufinal_tight.png']);
    
    ufigfinal.Position=fig_pos;
    set(gca,'FontSize',23);
    xlabel('x');
    ylabel('y');
    colorbar;
    set(gca,'XTick',0:(nx/numticks):nx);
    set(gca,'YTick',0:(nx/numticks):nx);
    set(gca,'XTickLabel',num2str((-L:2*L/numticks:L)'));
    set(gca,'YTickLabel',num2str((-L:2*L/numticks:L)'));
    title(['u, t=',num2str(T)],'FontSize',25);
    saveas(ufigfinal,[prefix,'_ufinal.png']);
    saveas(ufigfinal,[prefix,'_ufinal.fig']);
    
    vfigfinal=figure('Position',fig_pos_square,'color','w');
    imagesc(v, vrange);
    xlim([1,200]);
    ylim([1,200]);
    set(gca,'YDir','normal');
    axis image;
    colormap('hot');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(vfigfinal,[prefix,'_vfinal_tight.png']);
    
    vfigfinal.Position=fig_pos;
    set(gca,'FontSize',23);
    xlabel('x');
    ylabel('y');
    colorbar;
    set(gca,'XTick',0:(nx/numticks):nx);
    set(gca,'YTick',0:(nx/numticks):nx);
    set(gca,'XTickLabel',num2str((-L:2*L/numticks:L)'));
    set(gca,'YTickLabel',num2str((-L:2*L/numticks:L)'));
    title('v','FontSize',25);
    saveas(vfigfinal,[prefix,'_vfinal.png']);
    saveas(vfigfinal,[prefix,'_vfinal.fig']);
end
%% saving

if makegif
    ufinal = u;
    vfinal = v;
    save([prefix,'.mat'],'ufinal','vfinal','uyy_est','vyy_est', '-mat','-append');
end

end
