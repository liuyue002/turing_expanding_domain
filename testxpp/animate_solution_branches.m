function fig = animate_solution_branches(file,branchinds,n)
%Helper function to animate the PDE solution u as we trace along a branch
% branchinds: index of branches in the dat file
% n: which mode is this for
N=20;
nx=200;
fullbif=openfig('schnackenberg_fourier_auto_N=20.fig','invisible');
set(findall(gca, 'Type', 'Line'),'LineWidth',0.1);
databd = readbd(['schnackenberg_fourier_N=20_',num2str(file),'.dat']);

fig=figure('Position',[200,200,1200,500],'color','w');
giffile=['schnackenberg_fourier_N=20_mode=',num2str(n),'.gif'];
sfig1=subplot(1,2,1);
copyobj(allchild(get(fullbif, 'CurrentAxes')), sfig1);
hold on
for i=branchinds
    dat=databd.pts{1,i};
    type=databd.type{1,i};
    if strcmp(type,"ue")
        stability = '-k';
    elseif strcmp(type,"se")
        stability = '-r';
    end
    bifparam=dat(:,4);
    usum=sum(dat(:,7:2:7+N*2),2);
    plot(bifparam,usum,stability,'LineWidth',4);
end
xlim([0,50]);ylim([0,3]);
title(['Bifurcation diagram, n=',num2str(n),' branch']);
xlabel('L');
ylabel('u(0)');
dat=databd.pts{1,branchinds(1)};
for i=2:size(branchinds,2)
    dat=[dat;databd.pts{1,branchinds(i)}];
end
if n==1
    [~,I]=min(dat(:,4)); %I=478, start from the lower bif point
    dat=dat(I:I+1000,:); %cutoff, it just loops around
end
bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
numpts=size(dat,1);
skip=ceil(numpts/800);
marker=plot(bifparam(1),usum(1),'LineStyle','none','Marker','*','MarkerSize',20,'MarkerEdgeColor','b');
hold off

subplot(1,2,2);
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
for i =1:skip:size(dat,1)
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
end

