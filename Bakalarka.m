close all
clear

str = input("Do you want to change default setting? Y or enter: ",'s');
if ~isempty(str)
    newN = input("Set number of points: "); % check if integer
    newSmooth = input("Set smoothness of points: ");
    newRuns = input("Set number of runs: "); % don't have to set all
    newHeight = input("Set possible max height of points: ");
    % debug input
    [n, smoothVal, numRuns, height] = setVariables(newN, newSmooth, newRuns, newHeight);
    clear newN newSmooth newRuns; % was in local scope
else
    [n, smoothVal, numRuns, height] = setVariables(10, 0.2, 5, 50); % load from text file
end

[pointsArr, heightList] = generateEnvironment(smoothVal, n, height);
[r, c] = find(pointsArr > -1);
bools = pointsArr > -1;
coords = [r, c, pointsArr(bools)]; % almost x,y,z

triangArr = createTriangulations(r, c);
[lineList, lLwCoords] = createListOfLines(coords, triangArr);
sLWCoords = sortrows(lLwCoords, 'ascend');
smootherpoints = smoothenPlane(coords, sLWCoords); % EACH POINTS NEEDS TRIANGLE TO FIND HEIGHT
% CARTESIAN PRODUCT NOW
cartesianPoints = makeCartesianProduct(smootherpoints);
cartesianPoints = [cartesianPoints, ones(size(cartesianPoints,1), 1)];
[~, indx] = intersect(cartesianPoints(:,1:2), smootherpoints(:,1:2), 'rows');
cartesianPoints(indx, 3) = smootherpoints(:, 3);
disp(cartesianPoints);
% FIND BELONG TRIANGLE - for all triangles search remaining points
displayTriangulated(triangArr, heightList, r, c, cartesianPoints, pointsArr);
disp("DONE");

% Sets required variables and checks valid values
function [n, smoothVal, numRuns, height]  = setVariables(newN, newSmooth, newRuns, newHeight)
    n = checkValidity(newN, 4);
    smoothVal = checkValidity(newSmooth, 0.1);
    numRuns = checkValidity(newRuns, 1);
    height = checkValidity(newHeight, 0);
end

% Check input validity - NEED TO ADD STRINGS/FLOATS ...
function checkedValue = checkValidity(valueWanted, threshold)
    if valueWanted <= threshold
       checkedValue = threshold;
    else
       checkedValue = valueWanted;
    end
end

% Returns list of point heights with a 2D array with them marked
function [pointsArr, heightList] = generateEnvironment(smoothVal, n, height)
    % create appropriate 2d array
    nRC = 10 / (10 * smoothVal);
    pointsArr = zeros(nRC) - 1;
    nIdx = numel(pointsArr);
    noCornerIdx = [2:(nRC - 1),... % first row
        (nRC + 1):(nIdx - nRC),... % middle
        (nIdx - nRC + 2):(nIdx - 1)]; % last row
    cornerIdx = [1, nRC, nIdx - nRC + 1, nIdx];
    idx = sort([noCornerIdx(randperm(numel(noCornerIdx), n - 4)), cornerIdx]); % Matlab indexes from 1, maybe sort maybe no
    idx = [3,7,8,12,13,19, cornerIdx]; % testing crosses
    idx = [1, 3, nRC, 7, 8, 12, 13, 19, nIdx - nRC + 1, nIdx];
    heightRange = height + 1;
    heightList = randperm(heightRange, n) - 1;
    pointsArr(idx) = heightList; % can be 0
end

% Returns delaunay triangulation using delaunay library
function triangArr = createTriangulations(r, c)
    % create triangulation array from n and smoothVal
    triangArr = delaunay(r, c);
end

% Displays either 2D, 3D, or both
function displayTriangulated(triangArr, heightList, r, c, cartesianPoints, pointsArr)
    % display triangulated place
    str = input("Do you want 2D, 3D, or both visualizations? 2,3, or enter: ", 's');
    if (strcmp(str,"2"))
        triplot(triangArr, r, c);
        hold on
        plot(cartesianPoints(:, 1), cartesianPoints(:, 2), 'r.')
        hold off
    elseif (strcmp(str,"3"))    
        trimesh(triangArr, r, c, heightList);
        hold on
        plot(cartesianPoints(:, 1), cartesianPoints(:, 2), 'r.') % not enough
        hold off
    else
        triplot(triangArr, r, c);
        hold on
        plot(cartesianPoints(:, 1), cartesianPoints(:, 2), 'r.')
        hold off
        figure
        trimesh(triangArr, r, c, heightList); % indexes were manual and not ordered
        hold on
        plot3(cartesianPoints(:, 1), cartesianPoints(:, 2), cartesianPoints(:, 3), 'r.') % not enough
        hold off
    end
end

% Function returns lineList = list of lines with coords indices and
% lLwCoords = array of line coordinations (x1, y1, x2, y2)
function [lineList, lLwCoords] = createListOfLines(coords, triangArr)
    nTriangRows = size(triangArr, 1); % number triangulation lines
    maxLines = nTriangRows * 3; % number of triangles
    lineList = zeros(maxLines, 2);
    idxTaken = 1;
    for i = 1:nTriangRows
        noTriangles = idxTaken + 2;
        lineList(idxTaken:noTriangles, :) = nchoosek(triangArr(i,:), 2);
        idxTaken = noTriangles + 1;
    end
    lineList = unique(sort([lineList(:,1),lineList(:,2)], 2), 'rows');
    nOLines = size(lineList);
    lLwCoords = zeros(nOLines(1), 6); % preallc
    for i = 1:nOLines
        lLwCoords(i, :) = [coords(lineList(i,1), :), coords(lineList(i,2), :)]; %, coords(lineList(i,3), :)];
    end
end

function probDist = evaluateZeroSumGame(pointsUtils)
    % evaluate game
end

function smootherpoints = smoothenPlane(coords, sLWCoords)
    % create more points x, y % Gets called
    smootherpoints = coords
    for i = 1:size(coords, 1) % Runs correctly 10 times
        y =coords(i,:);
        x = ['Point ', num2str(y), ' has: '];
        disp(x);
        verticalInt = sLWCoords((sLWCoords(:,1)<coords(i,1) & sLWCoords(:,4)>coords(i,1)) | (sLWCoords(:,4)<coords(i,1) & sLWCoords(:,1)>coords(i,1)), :);
        %disp(verticalInt);
        % make points from his X same, calculate Y
        for j = 1:size(verticalInt, 1)
            if verticalInt(j,2) == verticalInt(j,5)
                newVerticalP = [coords(i,1), verticalInt(j,2), -1]
            else
                m = (verticalInt(j,2)-verticalInt(j,5)) / (verticalInt(j,1)-verticalInt(j,4));
                c = verticalInt(j,2) - m * verticalInt(j,1);
                newVerticalP = [coords(i,1), m * coords(i,1) + c, -1]
            end
            smootherpoints = cat(1, smootherpoints, newVerticalP);
        end
        horizontalInt = sLWCoords((sLWCoords(:,2)<coords(i,2) & sLWCoords(:,5)>coords(i,2)) | (sLWCoords(:,5)<coords(i,2) & sLWCoords(:,2)>coords(i,2)), :);
        %disp(horizontalInt);
        for j = 1:size(horizontalInt, 1)
            if horizontalInt(j,1) == horizontalInt(j,4)
                newHorizontalP = [horizontalInt(j,1), coords(i,2), -1]
            else
                m = (horizontalInt(j,2)-horizontalInt(j,5)) / (horizontalInt(j,1)-horizontalInt(j,4));
                c = horizontalInt(j,2) - m * horizontalInt(j,1);
                newHorizontalP = [(coords(i,2) - c) / m , coords(i,2), -1]
            end
            smootherpoints = cat(1, smootherpoints, newHorizontalP);
        end
        % make points from his X calculate, same Y
        %newHorizontalP =  [ , coords(i,2)];
        %smootherpoints = [smootherpoints, ...]
    end
    %disp(smootherpoints);
    smootherpoints = round(smootherpoints, 4);
    disp(smootherpoints); % OK
    smootherpoints = unique(smootherpoints, 'rows');
    disp(smootherpoints);
end

function cartesianPoints = makeCartesianProduct(smootherpoints)
    xn = unique(smootherpoints(:, 1), 'rows');
    yn = unique(smootherpoints(:, 2), 'rows');
    cartesianPoints = combvec(xn.', yn.').';
end

function newPoints = checkHorizontal()
end

function newPoints = checkVertical()
end
