function [fig] = plot_kymograph(uu, fig_pos,T,xrange,plot_title)
%Plot kymograph for PDE simulation results
% uu: matrix where uu(i,j) contains the value at t_i, x_j
% fig_pos: a vector of 4 numbers, specify location and size of figure
% T: total time of simulation
% xrange: range of x, i.e. [-L,L] or [0,L]
% plot_title: title for the figure
nFrame=size(uu,1);
nx=size(uu,2);
tTickNumber = 5;
tTick=(0:1/tTickNumber:1)*nFrame;
tTickLabel=num2cell((0:1/tTickNumber:1)*T);
format = cell(size(tTickLabel));
format(:)={'%.0f'};
tTickLabel=cellfun(@num2str,tTickLabel,format,'un',0);

xTickNumber = 10;
xTick=(0:1/xTickNumber:1)*nx;
xTickLabel=num2cell((0:1/xTickNumber:1)*(xrange(2)-xrange(1)) + xrange(1));
format = cell(size(xTickLabel));
format(:)={'%.0f'};
xTickLabel=cellfun(@num2str,xTickLabel,format,'un',0);

fig=figure('Position',fig_pos,'color','w');
axis([0 nFrame 0 nx]);
umax=max(max(uu(10:end,:)));
umax=ceil(umax*10)/10;
umin=min(min(uu(10:end,:)));
umin=floor(umin*10)/10;
ucolortick=[umin,umax];
imagesc(uu',ucolortick);
set(gca,'YDir','normal');
%colorbar('FontSize',40,'TickLabels',ucolortick,'Ticks',ucolortick);
colorbar('TickLabels',ucolortick,'Ticks',ucolortick);
set(gca,'XTick',tTick);
set(gca,'XTickLabel',tTickLabel);
set(gca,'YTick',xTick);
set(gca,'YTickLabel',xTickLabel);
xlabel('t');
ylabel('x');
title(plot_title);
%biggerFont(gca);
tightEdge(gca);

end

