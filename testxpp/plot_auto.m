addpath('./mdepitta');
N=20;
databd = readbd(['schnackenberg_fourier_N=',num2str(N),'.dat']);

%% plot all orbit
% dat=databd.pts{1,i} contains a single traced branch
% dat(:,1) is branch type, dat(:,2) is branch label
% dat(:,3) always 0?
% dat(:,4) is the bif param, dat(:,5) is secondary bif param
% dat(:,6) is period, dat(:,7) and on are the value of the ODE variables
% in the orders defined in the ODE file
% so here it's u0, v0, u1, v1... u[N], v[N]
% u[i] is column 7+2*i
% v[i] is column 8+2*i
% after that, it seems repeat u0,v0... u[N],v[N]
% after that, there's more. Maybe all the eigenvalues?
% databd.type{1,i}: looks like stability?
% ue=unstable eqlib, se = stable, uo= unstable orbit, so=stable orbit

fig=figure;
hold on
for i=1:size(databd.pts,2)
dat=databd.pts{1,i};
type=databd.type{1,i};
if strcmp(type,"ue")
    stability = '-k';
elseif strcmp(type,"se")
    stability = '-r';
else
    fprintf('Unknown type: %s\n',type);
    continue;
end

bifparam=dat(:,4);
usum=sum(dat(:,7:2:7+N*2),2);
line=plot(bifparam,usum,stability);
%line.Color = [line.Color, 0.2]; % make it more transparent
end
xlabel('L');
ylabel('u(0)');
xlim([0,50]);
ylim([0,3]);
title(['N=',num2str(N)]);
hold off
%saveas(fig,['schnackenberg_fourier_auto_N=',num2str(N),'.png']);