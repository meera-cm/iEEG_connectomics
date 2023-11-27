% Plot Correlation Matrix

% Produce the input lower triangular matrix data
C = -1 + 2.*rand(12,12);
C = tril(C,-1); % tril returns the lower triangular portion of the data
C(logical(eye(size(C)))) = 1;
% Set [min,max] value of C to scale colors
clrLim = [-1,1];
% load('CorrColormap.mat') % Uncomment for custom CorrColormap
% Set the  [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.1, 1];
myLabel = {'ICA','Elev','Pr','Rmax','Rmin','Srad','Wspd','Tmin','Tmax','VPD','ET_o','AW'};
% Compute center of each circle
% This assumes the x and y values were not entered in imagesc()
x = 1 : 1 : size(C,2); % x edges
y = 1 : 1 : size(C,1); % y edges
[xAll, yAll] = meshgrid(x,y);
xAll(C==0)=nan; % eliminate cordinates for zero correlations
% Set color of each rectangle
% Set color scale
cmap = jet(256);
% cmap = CorrColormap; % Uncomment for CorrColormap
Cscaled = (C - clrLim(1))/range(clrLim); % always [0:1]
colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size between [0 1]
Cscaled = (abs(C) - 0)/1;
diamSize = Cscaled * range(diamLim) + diamLim(1);
% Create figure
fh = figure();
ax = axes(fh);
hold(ax,'on')
colormap(ax,'jet');
% colormap(CorrColormap) %Uncomment for CorrColormap
tickvalues = 1:length(C);
x = zeros(size(tickvalues));
text(x, tickvalues, myLabel, 'HorizontalAlignment', 'right');
x(:) = length(C)+1;
text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',90);
% Create circles
theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:),'LineStyle','none'),1:numel(xAll));
axis(ax,'equal')
axis(ax,'tight')
set(ax,'YDir','Reverse')
colorbar()
caxis(clrLim);
axis off
