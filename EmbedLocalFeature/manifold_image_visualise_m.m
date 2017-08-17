function manifold_image_visualize_m(X, Y, visualiseFunction, axesWidth, varargin)

% visualized imeage embedding in low dimensional(2D) manifold
% modified from gplvmStaticImageVisualise 
% by Chan-Su Lee Aug. 28, 2005
%
% GPLVMSTATICIMAGEVISUALISE Generate a scatter plot of the images without overlap.
%
% gplvmStaticImageVisualise(model, visualiseFunction, axesWidth, varargin)

% Copyright (c) 2004 Neil D. Lawrence
% File version 1.3, Thu Jul  8 15:31:53 2004
% GPLVM toolbox version 2.011



% set random seeds
randn('seed', 1e5)
rand('seed', 1e5)

colordef white 

% % Turn Y into grayscale
% h = figure
% clf
% % Create the plot for the data
% clf

% % plot(X(:,1), X(:,2), 'ro-')
% plot(X(:,1), X(:,2), 'ro-')
% axis off

% posit = [0.05 0.05 0.9 0.9];
% posit = [0.0 0.0 1.0 1.0];
% plot(X(:,1), X(:,2), 'ro-')
% axis on

% ax = axes('position', posit);
% 
% xLim = [min(X(:, 1)) max(X(:, 1))];
% yLim = [min(X(:, 2)) max(X(:, 2))];
% set(ax, 'xLim', xLim);
% set(ax, 'yLim', yLim);
% 
% set(ax, 'fontname', 'arial');
% set(ax, 'fontsize', 20);

%     plot(X(:,1), X(:,2), 'ro-')
%     pi_str = int2str(pindex(i));
%     K_str = int2str(K);
%     title(['Subject = ' pi_str ' K = ' K_str ' D=3'])
    
    
[plotAxes, data] = manifold_scatter_plot(X, []);

xLim = get(plotAxes, 'xLim');
yLim = get(plotAxes, 'yLim');
posit = get(plotAxes, 'position');

widthVal = axesWidth*(xLim(2) - xLim(1))/posit(3);
heightVal = axesWidth*(yLim(2) - yLim(1))/posit(4);
numData = size(X, 1);




visitOrder = randperm(numData);
initVisitOrder = visitOrder;

% newPointsX = [];
% for j=1:numData 
%     point = invGetNormAxesPoint(X(j, :), ax);
%     newPointsX = [newPointsX; point];
% end
% 
% plot(newPointsX(:,1), newPointsX(:,2), 'ro-')

%     digitAxes(i) =  axes('position', ...
% 			 [x - axesWidth/2 ...
% 		    y - axesWidth/2 ...
% 		    axesWidth ...
% 		    axesWidth]);
%     handle = feval(visualiseFunction, Y(i, :), varargin{:});
%     colormap gray
%     axis image
%     axis off
%     
%     removePoints = find(abs(X(visitOrder, 1) - X(i, 1)) < widthVal/2 ...
% 			&  abs(X(visitOrder, 2) - X(i, 2)) < heightVal);
%     visitOrder(removePoints) = [];    
% end


% Plot the images
while ~isempty(visitOrder)
  i = visitOrder(1);
  if X(i, 1) > xLim(1) & X(i, 1) < xLim(2) ...
    & X(i, 2) > yLim(1) & X(i, 2) < yLim(2)
%     point = invGetNormAxesPoint(X(i, :), ax);
    point = invGetNormAxesPoint(X(i, :), plotAxes);
    x = point(1);
    y = point(2);
    
    digitAxes(i) =  axes('position', ...
			 [x - axesWidth/2 ...
		    y - axesWidth/2 ...
		    axesWidth ...
		    axesWidth]);
    handle = feval(visualiseFunction, Y(i, :), varargin{:});
    colormap gray
    axis image
    axis off
    
    %removePoints = find(abs(X(visitOrder, 1) - X(i, 1)) < widthVal/2 ...
	%		&  abs(X(visitOrder, 2) - X(i, 2)) < heightVal);
    %visitOrder(removePoints) = [];    
    visitOrder(1) = [];
    
  else
    visitOrder(1) = [];
  end
end
    
% plot(X(:,1), X(:,2), 'ro-')
    
% set(ax, 'xlim', xLim);
% set(ax, 'ylim', yLim);
% set(ax, 'visible', 'off');
set(plotAxes, 'xlim', xLim);
set(plotAxes, 'ylim', yLim);
set(plotAxes, 'visible', 'off');
%ticks = [-4 -2 0 2 4];
%set(plotAxes, 'xtick', ticks)
%set(plotAxes, 'ytick', ticks)
