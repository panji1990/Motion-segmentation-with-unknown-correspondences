function [ax, data] = manifold_scatter_plot(X, YLbls)
% function [ax, data] = manifold_scatter_plot(X, YLbls)
% modified from gplvmScatterPlot(model, YLbls);
% by Chan-Su Lee  Aug. 29, 2005
%
% GPLVMSCATTERPLOT 2-D scatter plot of the latent points.
%
% [ax, data] = gplvmScatterPlot(model, YLbls);

% Copyright (c) 2004 Neil D. Lawrence
% File version 1.1, Mon Jul 19 08:40:00 2004
% GPLVM toolbox version 2.011



if isempty(YLbls)
  symbol = [];
else
  symbol = getSymbols(size(YLbls,2));
end

% if nargin<3
%     title_str ='Manifold plot with images';
% end

% x1 = linspace(min(model.X(:, 1))*1.1, max(model.X(:, 1))*1.1, 30);
% x2 = linspace(min(model.X(:, 2))*1.1, max(model.X(:, 2))*1.1, 30);
% [X1, X2] = meshgrid(x1, x2);
% XTest = [X1(:), X2(:)];
% [mu, varsigma] = ivmPosteriorMeanVar(model, XTest);
%   

figure
clf
% Create the plot for the data
clf
% title(title_str)
ax = axes('position', [0.05 0.05 0.9 0.9]);
% ax = axes('position', [0.0 0.0 1.0 1.0]);
hold on
% [c, h] = contourf(X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 128); 
% shading flat
% colormap gray;
%colorbar
border_offset = 0.2;
data = gplvmtwoDPlot(X, YLbls, symbol);

xLim = [min(X(:, 1)-border_offset) max(X(:, 1))+border_offset];
yLim = [min(X(:, 2))-border_offset max(X(:, 2))+border_offset];
set(ax, 'xLim', xLim);
set(ax, 'yLim', yLim);

set(ax, 'fontname', 'arial');
set(ax, 'fontsize', 20);

function returnVal = gplvmtwoDPlot(X, label, symbol)

% GPLVMTWODPLOT Helper function for plotting the labels in 2-D.

returnVal = [];

if ~isempty(label)
  for i = 1:size(X, 1)
    labelNo = find(label(i, :));
    try 
      returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo},'LineWidth',2)];
    catch
      if strcmp(lasterr, 'Index exceeds matrix dimensions.')
	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
      end
    end
  end
else
  %returnVal = plot(X(:, 1), X(:, 2), 'rx-','LineWidth',2);
  returnVal = plot(X(:, 1), X(:, 2), 'r.','LineWidth',2);
end