addpath('./mdepitta');
N=20;
nx=200;

%% mode n=1 exploration, branch is plot _1_7
file=1;
branch=7;
databd = readbd(['schnackenberg_fourier_N=20_',num2str(file),'.dat']);
type=databd.type{1,branch};
dat=databd.pts{1,branch};
[~,I]=min(dat(:,4)); %I=478, start from the lower bif point
dat=dat(I:I+1000,:); %cutoff, it just loops around

figure;
plot(dat(:,4)); % how the bif param goes
xlabel('Index');
ylabel('L');


%% mode n=1 plotting
fig=figure('Position',[200,200,1200,500],'color','w');
giffile='schnackenberg_fourier_N=20_mode=1.gif';
sfig1=subplot(1,2,1);

if strcmp(type,"ue")
    stability = '-k';
elseif strcmp(type,"se")
    stability = '-r';
end
bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
hold on
plot(bifparam,usum,stability);
marker=plot(bifparam(1),usum(1),'LineStyle','none','Marker','*','MarkerSize',10,'MarkerEdgeColor','b');
hold off
xlim([0,50]);ylim([0,3]);
title('Bifurcation diagram, N=1 branch');
xlabel('L');
ylabel('u(0)');

sfig2=subplot(1,2,2);
L=dat(1,4);
coefs=dat(1,7:2:7+N*2);
x=linspace(0,L,nx);
u=zeros(1,nx);
for n=0:N
    u=u+coefs(n+1)*cos((n*pi/L)*x);
end
uplot=plot(x,u);
ylim([0,3]);
plottitle=title(['The solution, L=',num2str(L,'%.3f')]);
xlabel('x');
ylabel('u');
[imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
for i =1:2:size(dat,1)
    L=dat(i,4);
    coefs=dat(i,7:2:7+N*2);
    x=linspace(0,L,nx);
    u=zeros(1,nx);
    for n=0:N
        u=u+coefs(n+1)*cos((n*pi/L)*x);
    end
    uplot.XData=x;
    uplot.YData=u;
    marker.XData=L;
    marker.YData=usum(i);
    xlim([0,L]);
    plottitle.String=['The solution, L=',num2str(L,'%.3f')];
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
    imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
end

%% mode n=2 plotting, branch is plot _1_23 through 26
file=1;
databd = readbd(['schnackenberg_fourier_N=20_',num2str(file),'.dat']);
fig=figure('Position',[200,200,1200,500],'color','w');
giffile='schnackenberg_fourier_N=20_mode=2.gif';
sfig1=subplot(1,2,1);
hold on
for i=23:26
    dat=databd.pts{1,i};
    type=databd.type{1,i};
    if strcmp(type,"ue")
        stability = '-k';
    elseif strcmp(type,"se")
        stability = '-r';
    end
    bifparam=dat(:,4);
    usum=sum(dat(:,7:2:7+N*2),2);
    plot(bifparam,usum,stability);
end
xlim([0,50]);ylim([0,3]);
title('Bifurcation diagram, N=2 branch');
xlabel('L');
ylabel('u(0)');
dat=[databd.pts{1,23};databd.pts{1,24};databd.pts{1,25};databd.pts{1,26}];
bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
marker=plot(bifparam(1),usum(1),'LineStyle','none','Marker','*','MarkerSize',10,'MarkerEdgeColor','b');
hold off

sfig2=subplot(1,2,2);
L=dat(1,4);
coefs=dat(1,7:2:7+N*2);
x=linspace(0,L,nx);
u=zeros(1,nx);
for n=0:N
    u=u+coefs(n+1)*cos((n*pi/L)*x);
end
uplot=plot(x,u);
ylim([0,3]);
plottitle=title(['The solution, L=',num2str(L,'%.3f')]);
xlabel('x');
ylabel('u');
[imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
for i =1:size(dat,1)
    L=dat(i,4);
    coefs=dat(i,7:2:7+N*2);
    x=linspace(0,L,nx);
    u=zeros(1,nx);
    for n=0:N
        u=u+coefs(n+1)*cos((n*pi/L)*x);
    end
    uplot.XData=x;
    uplot.YData=u;
    marker.XData=L;
    marker.YData=usum(i);
    xlim([0,L]);
    plottitle.String=['The solution, L=',num2str(L,'%.3f')];
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
    imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
end

%% mode n=3 plotting, branch is plot _1_42-46

databd = readbd('schnackenberg_fourier_N=20_1.dat');
fig=figure('Position',[200,200,1200,500],'color','w');
giffile='schnackenberg_fourier_N=20_mode=3.gif';
sfig1=subplot(1,2,1);
hold on
for i=42:46
    dat=databd.pts{1,i};
    type=databd.type{1,i};
    if strcmp(type,"ue")
        stability = '-k';
    elseif strcmp(type,"se")
        stability = '-r';
    end
    bifparam=dat(:,4);
    usum=sum(dat(:,7:2:7+N*2),2);
    plot(bifparam,usum,stability);
end
xlim([0,50]);ylim([0,3]);
title('Bifurcation diagram, n=3 branch');
xlabel('L');
ylabel('u(0)');
dat=[databd.pts{1,42};databd.pts{1,43};databd.pts{1,44};databd.pts{1,45};databd.pts{1,46}];
bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
marker=plot(bifparam(1),usum(1),'LineStyle','none','Marker','*','MarkerSize',10,'MarkerEdgeColor','b');
hold off

sfig2=subplot(1,2,2);
L=dat(1,4);
coefs=dat(1,7:2:7+N*2);
x=linspace(0,L,nx);
u=zeros(1,nx);
for n=0:N
    u=u+coefs(n+1)*cos((n*pi/L)*x);
end
uplot=plot(x,u);
ylim([0,3]);
plottitle=title(['The solution, L=',num2str(L,'%.3f')]);
xlabel('x');
ylabel('u');
[imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
for i =1:3:size(dat,1)
    L=dat(i,4);
    coefs=dat(i,7:2:7+N*2);
    x=linspace(0,L,nx);
    u=zeros(1,nx);
    for n=0:N
        u=u+coefs(n+1)*cos((n*pi/L)*x);
    end
    uplot.XData=x;
    uplot.YData=u;
    marker.XData=L;
    marker.YData=usum(i);
    xlim([0,L]);
    plottitle.String=['The solution, L=',num2str(L,'%.3f')];
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
    imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
end

%% mode n=4 plotting, branch is plot _1_, 56-60

databd = readbd('schnackenberg_fourier_N=20_1.dat');
fig=figure('Position',[200,200,1200,500],'color','w');
giffile='schnackenberg_fourier_N=20_mode=4.gif';
sfig1=subplot(1,2,1);
hold on
for i=56:60
    dat=databd.pts{1,i};
    type=databd.type{1,i};
    if strcmp(type,"ue")
        stability = '-k';
    elseif strcmp(type,"se")
        stability = '-r';
    end
    bifparam=dat(:,4);
    usum=sum(dat(:,7:2:7+N*2),2);
    plot(bifparam,usum,stability);
end
xlim([0,50]);ylim([0,3]);
title('Bifurcation diagram, n=4 branch');
xlabel('L');
ylabel('u(0)');
dat=[databd.pts{1,56};databd.pts{1,57};databd.pts{1,58};databd.pts{1,59};databd.pts{1,60}];
bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
marker=plot(bifparam(1),usum(1),'LineStyle','none','Marker','*','MarkerSize',10,'MarkerEdgeColor','b');
hold off

sfig2=subplot(1,2,2);
L=dat(1,4);
coefs=dat(1,7:2:7+N*2);
x=linspace(0,L,nx);
u=zeros(1,nx);
for n=0:N
    u=u+coefs(n+1)*cos((n*pi/L)*x);
end
uplot=plot(x,u);
ylim([0,3]);
plottitle=title(['The solution, L=',num2str(L,'%.3f')]);
xlabel('x');
ylabel('u');
[imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
for i =1:3:size(dat,1)
    L=dat(i,4);
    coefs=dat(i,7:2:7+N*2);
    x=linspace(0,L,nx);
    u=zeros(1,nx);
    for n=0:N
        u=u+coefs(n+1)*cos((n*pi/L)*x);
    end
    uplot.XData=x;
    uplot.YData=u;
    marker.XData=L;
    marker.YData=usum(i);
    xlim([0,L]);
    plottitle.String=['The solution, L=',num2str(L,'%.3f')];
    [imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
    imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
end

%% mode n=5 plotting, branch is plot _2_, 2-6

databd = readbd('schnackenberg_fourier_N=20_2.dat');
fig=figure('Position',[200,200,1200,500],'color','w');
giffile='schnackenberg_fourier_N=20_mode=5.gif';
sfig1=subplot(1,2,1);
hold on
for i=2:6
    dat=databd.pts{1,i};
    type=databd.type{1,i};
    if strcmp(type,"ue")
        stability = '-k';
    elseif strcmp(type,"se")
        stability = '-r';
    end
    bifparam=dat(:,4);
    usum=sum(dat(:,7:2:7+N*2),2);
    plot(bifparam,usum,stability);
end
xlim([0,50]);ylim([0,3]);
title('Bifurcation diagram, n=5 branch');
xlabel('L');
ylabel('u(0)');
dat=[databd.pts{1,2};databd.pts{1,3};databd.pts{1,4};databd.pts{1,5};databd.pts{1,6}];
bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
marker=plot(bifparam(1),usum(1),'LineStyle','none','Marker','*','MarkerSize',10,'MarkerEdgeColor','b');
hold off

sfig2=subplot(1,2,2);
L=dat(1,4);
coefs=dat(1,7:2:7+N*2);
x=linspace(0,L,nx);
u=zeros(1,nx);
for n=0:N
    u=u+coefs(n+1)*cos((n*pi/L)*x);
end
uplot=plot(x,u);
ylim([0,3]);
plottitle=title(['The solution, L=',num2str(L,'%.3f')]);
xlabel('x');
ylabel('u');
[imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
for i =1:5:size(dat,1)
    L=dat(i,4);
    coefs=dat(i,7:2:7+N*2);
    x=linspace(0,L,nx);
    u=zeros(1,nx);
    for n=0:N
        u=u+coefs(n+1)*cos((n*pi/L)*x);
    end
    uplot.XData=x;
    uplot.YData=u;
    marker.XData=L;
    marker.YData=usum(i);
    xlim([0,L]);
    plottitle.String=['The solution, L=',num2str(L,'%.3f')];
    imind = rgb2ind(frame2im(getframe(fig)),cm);
    imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
end

%% mode n=6 plotting, branch is plot _2_, 24-28

databd = readbd('schnackenberg_fourier_N=20_2.dat');
fig=figure('Position',[200,200,1200,500],'color','w');
giffile='schnackenberg_fourier_N=20_mode=6.gif';
sfig1=subplot(1,2,1);
hold on
for i=24:28
    dat=databd.pts{1,i};
    type=databd.type{1,i};
    if strcmp(type,"ue")
        stability = '-k';
    elseif strcmp(type,"se")
        stability = '-r';
    end
    bifparam=dat(:,4);
    usum=sum(dat(:,7:2:7+N*2),2);
    plot(bifparam,usum,stability);
end
xlim([0,50]);ylim([0,3]);
title('Bifurcation diagram, n=6 branch');
xlabel('L');
ylabel('u(0)');
dat=[databd.pts{1,24};databd.pts{1,25};databd.pts{1,26};databd.pts{1,27};databd.pts{1,28}];
bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
marker=plot(bifparam(1),usum(1),'LineStyle','none','Marker','*','MarkerSize',10,'MarkerEdgeColor','b');
hold off

sfig2=subplot(1,2,2);
L=dat(1,4);
coefs=dat(1,7:2:7+N*2);
x=linspace(0,L,nx);
u=zeros(1,nx);
for n=0:N
    u=u+coefs(n+1)*cos((n*pi/L)*x);
end
uplot=plot(x,u);
ylim([0,3]);
plottitle=title(['The solution, L=',num2str(L,'%.3f')]);
xlabel('x');
ylabel('u');
[imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
imwrite(imind,cm,giffile,'gif', 'Loopcount',inf);
for i =1:5:size(dat,1)
    L=dat(i,4);
    coefs=dat(i,7:2:7+N*2);
    x=linspace(0,L,nx);
    u=zeros(1,nx);
    for n=0:N
        u=u+coefs(n+1)*cos((n*pi/L)*x);
    end
    uplot.XData=x;
    uplot.YData=u;
    marker.XData=L;
    marker.YData=usum(i);
    xlim([0,L]);
    plottitle.String=['The solution, L=',num2str(L,'%.3f')];
    imind = rgb2ind(frame2im(getframe(fig)),cm);
    imwrite(imind,cm,giffile,'gif','WriteMode','append','DelayTime',0);
end