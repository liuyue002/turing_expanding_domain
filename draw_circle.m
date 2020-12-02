function circle_object = draw_circle(xc, yc, r)
%Draw a circle on the current figure
hold on
theta = linspace(0, 2*pi, 200);
circle_object = plot(xc + r*cos(theta), yc+r*sin(theta), 'b-');
hold off
end

