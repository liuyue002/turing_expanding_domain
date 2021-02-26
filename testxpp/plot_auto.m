addpath('./mdepitta');
N=20;
databd = readbd(['schnackenberg_fourier_N=',num2str(N),'.dat']);
unorm = "L2"; % "boundary" or "L2" or "W12"

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
for datind=1:size(databd.pts,2)
    dat=databd.pts{1,datind};
    type=databd.type{1,datind};
    if strcmp(type,"ue")
        stability = '-k';
    elseif strcmp(type,"se")
        stability = '-r';
    else
        fprintf('Unknown type: %s\n',type);
        continue;
    end
    
    bifparam=dat(:,4);
    if unorm=="boundary"
        usum=sum(dat(:,7:2:7+N*2),2); % this is u(0)
        line=plot(bifparam,usum,stability);
    elseif unorm=="L2"
        L2norm=zeros(size(dat,1),1);
        nx=200;
        for i=1:size(dat,1)
            L=bifparam(i);
            coefs=dat(i,7:2:7+N*2);
            x=linspace(0,L,nx);
            u=zeros(1,nx);
            for n=0:N
                u=u+coefs(n+1)*cos((n*pi/L)*x);
            end
            L2norm(i)=sqrt(sum(u.^2)/nx);
        end
        line=plot(bifparam,L2norm,stability);
    else
        %unorm=="W12"
        W12norm=zeros(size(dat,1),1);
        nx=200;
        for i=1:size(dat,1)
            L=bifparam(i);
            coefs=dat(i,7:2:7+N*2);
            x=linspace(0,L,nx);
            u=zeros(1,nx);
            uder=zeros(1,nx);
            for n=0:N
                u=u+coefs(n+1)*cos((n*pi/L)*x);
                uder=uder-coefs(n+1)*(n*pi/L)*sin((n*pi/L)*x);
            end
            W12norm(i)=sqrt((sum(u.^2)+sum(uder.^2))/nx);
        end
        line=plot(bifparam,W12norm,stability);
    end
    
    %line.Color = [line.Color, 0.2]; % make it more transparent
end

xlabel('L');
if unorm=="boundary"
    ylabel('u(0)');
    xlim([0,50]);
    ylim([0,3]);
elseif unorm=="L2"
    ylabel('||u||_2 /L');
    xlim([0,50]);
    ylim([1.4,1.7]);
else
    %unorm=="W12"
    ylabel('||u||_{1,2} /L');
    xlim([0,50]);
    ylim([1.4,1.7]);
end

title(['N=',num2str(N)]);
hold off
saveas(fig,strcat('schnackenberg_fourier_auto_N=',num2str(N),'_unorm=',unorm,'.png'));

