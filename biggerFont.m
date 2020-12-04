function [ ] = biggerFont( ax )
%Makes all text in figure bigger and lines thicker for publication
drawnow;
set(ax,'FontSize', 40);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
end
