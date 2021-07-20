function [fig] = plot_kymograph(uu, fig_pos,T,xrange,urange,plot_title,forprint)
%Plot kymograph for PDE simulation results
% uu: matrix where uu(i,j) contains the value at t_i, x_j
% fig_pos: a vector of 4 numbers, specify location and size of figure
% T: total time of simulation
% xrange: range of x, i.e. [-L,L] or [0,L]
% urange: range for color. Leave at NaN to pick the range automatically
% plot_title: title for the figure
% forprint: whether the figure is for publication (larger fonts, etc)
nFrame=size(uu,1);
nx=size(uu,2);

if ~forprint
    tTickNumber = 5;
    xTickNumber = 10;
else
    tTickNumber = 2;
    xTickNumber = 2;
end

tTick=(0:1/tTickNumber:1)*nFrame;
tTick(1)=1;
tTickLabel=num2cell((0:1/tTickNumber:1)*T);
format = cell(size(tTickLabel));
format(:)={'%.0f'};
tTickLabel=cellfun(@num2str,tTickLabel,format,'un',0);

xTick=(0:1/xTickNumber:1)*nx;
xTick(1)=1;
xTickLabel=num2cell((0:1/xTickNumber:1)*(xrange(2)-xrange(1)) + xrange(1));
format = cell(size(xTickLabel));
format(:)={'%.0f'};
xTickLabel=cellfun(@num2str,xTickLabel,format,'un',0);

fig=figure('Position',fig_pos,'color','w');
axis([0 nFrame 0 nx]);
if isnan(urange)
    umax=max(max(uu(10:end,:)));
    umax=ceil(umax*10)/10;
    umin=min(min(uu(10:end,:)));
    umin=floor(umin*10)/10;
    urange=[umin,umax];
end
imagesc(uu',urange);
set(gca,'YDir','normal');
%colorbar('FontSize',40,'TickLabels',ucolortick,'Ticks',ucolortick);
set(gca,'XTick',tTick);
set(gca,'XTickLabel',tTickLabel);
set(gca,'YTick',xTick);
set(gca,'YTickLabel',xTickLabel);
xlabel('t');
ylabel('x');
title(plot_title);
if ~forprint
    colorbar('TickLabels',urange,'Ticks',urange);
else
    biggerFont(gca);
end
tightEdge(gca);

end

