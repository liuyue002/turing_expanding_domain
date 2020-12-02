function update_circle(circle_object, xc, yc, r)
%Update a circle object on a plot
theta = linspace(0, 2*pi, 200);
circle_object.XData=xc + r*cos(theta);
circle_object.YData=yc+r*sin(theta);
end

