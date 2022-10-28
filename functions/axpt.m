function axesPosition = axpt(nX, nY, xPos, yPos, positionVector, interval)
%axpt gives position vector for the function axes
%
%   nX: number of column
%   nY: number of row
%   xPos: x position of subplot. can be scalar or row vector (ex. 4 or [3:5])
%       the most left column is 1. the most right column is nX
%   yPos: y position of subplot. can be scalar or row vector
%       the most upper row is 1. the most lower row is nY
%   interval: interval between col and row. it is proportional
%   to width and height of figure. should be two row vector ([0.05 0.05])
%   positionVector: default [0.1 0.1 0.85 0.85];


narginchk(4,6);
if nargin <= 5
    interval = [0.05 0.05];
    if nargin <= 4
        positionVector = [0.1 0.1 0.85 0.85];
    end  
end
if isempty(positionVector)
    positionVector = [0.1 0.1 0.85 0.85];
end

xInterval = interval(1)*positionVector(3);
yInterval = interval(2)*positionVector(4);
dX = (positionVector(3) - (nX-1)*xInterval) / nX;
dY = (positionVector(4) - (nY-1)*yInterval) / nY;

xMin = positionVector(1) + (min(xPos)-1)*(dX+xInterval);
xMax = positionVector(1) + (max(xPos)-1)*(dX+xInterval) + dX;

yMin = positionVector(2) + (nY-max(yPos))*(dY+yInterval);
yMax = positionVector(2) + (nY-min(yPos))*(dY+yInterval) + dY;

axesPosition = [xMin, yMin, xMax-xMin, yMax-yMin];