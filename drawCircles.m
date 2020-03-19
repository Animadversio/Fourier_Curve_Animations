function drawCircles(centers, radii, varargin)
thetaResolution = 2; 
theta=(0:thetaResolution:360)*pi/180;
x = radii*cos(theta)+centers(:,1);
x = cat(1,x',nan(1,length(radii)));
x = x(:);
y = radii*sin(theta)+centers(:,2);
y = cat(1,y',nan(1,length(radii)));
y = y(:);
line(x,y,'Color',[0.3010, 0.7450, 0.9330],"LineWidth",0.5,'LineStyle','-.')%[0, 0.2235, 0.3705]
% ...
%    'Color',options.Color, ...
%    'LineWidth',options.LineWidth, ...
%    'LineStyle',options.LineStyle);
end
