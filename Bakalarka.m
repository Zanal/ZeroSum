close all
clear

%runThesis();

% Function for 
%function runThesis()
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
coords = [r, c, pointsArr(bools)]; % x,y,z

triangArr = createTriangulations(r, c);
[lineList, lLwCoords] = createListOfLines(coords, triangArr);
sLWCoords = sortrows(lLwCoords, 'ascend');
smootherpoints = smoothenPlane(coords, sLWCoords); % EACH POINTS NEEDS TRIANGLE TO FIND HEIGHT
cartesianPoints = makeCartesianProduct(smootherpoints);
cartesianPoints = [cartesianPoints, -ones(size(cartesianPoints,1), 1)];
[~, indx] = intersect(cartesianPoints(:,1:2), smootherpoints(:,1:2), 'rows');
cartesianPoints(indx, 3) = smootherpoints(:, 3);
cartesianPoints = calculateCPHeights(triangArr, coords, cartesianPoints);
disp(cartesianPoints);
% FIND BELONG TRIANGLE - for all triangles search remaining points + height
displayTriangulated(triangArr, heightList, r, c, cartesianPoints);
disp("DONE");

%end

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
    idx = [1, 3, nRC, 7, 8, 12, 13, 19, nIdx - nRC + 1, nIdx]; % TESTING
    heightRange = height + 1;
    heightList = randperm(heightRange, n) - 1;
    heightList = [4,41,34,15,44,1,19,16,32,33]; % TESTING
    pointsArr(idx) = heightList; % can be 0
end

% Returns delaunay triangulation using delaunay library
function triangArr = createTriangulations(r, c)
    % create triangulation array from n and smoothVal
    triangArr = delaunay(r, c);
end

% Displays either 2D, 3D, or both
function displayTriangulated(triangArr, heightList, r, c, cartesianPoints)
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
        plot3(cartesianPoints(:, 1), cartesianPoints(:, 2), cartesianPoints(:, 3), 'r.') % not enough
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

% Find new intersections
function smootherpoints = smoothenPlane(coords, sLWCoords)
    % create more points x, y % Gets called
    smootherpoints = coords;
    for i = 1:size(coords, 1) % Runs correctly 10 times
        verticalInt = sLWCoords((sLWCoords(:,1)<coords(i,1) & sLWCoords(:,4)>coords(i,1)) | ...
                        (sLWCoords(:,4)<coords(i,1) & sLWCoords(:,1)>coords(i,1)), :);
        % make points from his X same, calculate Y
        for j = 1:size(verticalInt, 1)
            if verticalInt(j,2) == verticalInt(j,5)
                newVerticalP = [coords(i,1), verticalInt(j,2), -1];
            else
                m = (verticalInt(j,2)-verticalInt(j,5)) / (verticalInt(j,1)-verticalInt(j,4));
                c = verticalInt(j,2) - m * verticalInt(j,1);
                newVerticalP = [coords(i,1), m * coords(i,1) + c, -1];
            end
            newVerticalP(3) = calculateHeights(verticalInt(j,:), newVerticalP);
            smootherpoints = cat(1, smootherpoints, newVerticalP);
        end
        horizontalInt = sLWCoords((sLWCoords(:,2)<coords(i,2) & sLWCoords(:,5)>coords(i,2)) | ...
                            (sLWCoords(:,5)<coords(i,2) & sLWCoords(:,2)>coords(i,2)), :);
        for j = 1:size(horizontalInt, 1)
            if horizontalInt(j,1) == horizontalInt(j,4)
                newHorizontalP = [horizontalInt(j,1), coords(i,2), -1];
            else
                m = (horizontalInt(j,2)-horizontalInt(j,5)) / (horizontalInt(j,1)-horizontalInt(j,4));
                c = horizontalInt(j,2) - m * horizontalInt(j,1);
                newHorizontalP = [(coords(i,2) - c) / m , coords(i,2), -1];
            end
            newHorizontalP(3) = calculateHeights(horizontalInt(j,:), newHorizontalP);
            smootherpoints = cat(1, smootherpoints, newHorizontalP);
        end
    end
    smootherpoints = round(smootherpoints, 4);
    smootherpoints = unique(smootherpoints, 'rows');
end

% Function used to calculate height of points found by intersection
function calcedHeight = calculateHeights(linePoints, foundPoint)
    if linePoints(3) < linePoints(6)
        smallerCoords = linePoints(1:2);
        smallerH = linePoints(3);
    else
        smallerCoords = linePoints(4:5);
        smallerH = linePoints(6);
    end
    calcedHeight = smallerH + norm(foundPoint(1:2)-smallerCoords) / ...
        norm(linePoints(1:2)-linePoints(4:5)) * (abs(linePoints(3)-linePoints(6)));
end

% Creates Cartesian product from unique x and y coordinates
function cartesianPoints = makeCartesianProduct(smootherpoints)
    xn = unique(smootherpoints(:, 1), 'rows');
    yn = unique(smootherpoints(:, 2), 'rows');
    cartesianPoints = combvec(xn.', yn.').';
end

function cartesianPoints = calculateCPHeights(triangArr, coords, cartesianPoints)
    % for all pointss with -1 (similarly like putting height in runner)
    % go through triangles to see if they belong such as c1,c2,c3 in py
    % calculate height using linspace? like py & assign
    noHeightIdx = find(cartesianPoints(:,3) == -1)
    for i = 1:length(noHeightIdx)
       % triangArr - idx of points in coords
       for j = 1:length(triangArr)
            %cartesianPoints(j, 1:2) % x,y
            p1 = coords(triangArr(j,1), :); % x, y 1
            p2 = coords(triangArr(j,2), :); % x,y 2
            p3 = coords(triangArr(j,3), :); % x,y 3
            chp = cartesianPoints(noHeightIdx(i), 1:2); % x,y checked
            c1 = (p2(1) - p1(1)) * (chp(2) - p1(2)) - (p2(2) - p1(2)) * (chp(1) - p1(1));
            c2 = (p3(1) - p2(1)) * (chp(2) - p2(2)) - (p3(2) - p2(2)) * (chp(1) - p2(1));
            c3 = (p1(1) - p3(1)) * (chp(2) - p3(2)) - (p1(2) - p3(2)) * (chp(1) - p3(1)); % OK
            disp(j);
            if (c1 <= 0 && c2 <= 0 && c3 <= 0) || (c1 >= 0 && c2 >= 0 && c3 >= 0)
                dd1 = ((p2(1)-p1(1))*(p3(3)-p1(3)) - (p3(1)-p1(1))*(p2(3)-p1(3))); % function for all triangles to avoid repetative calculation
                ds1 = ((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)));
                dd2 = ((p2(2)-p1(2))*(p3(3)-p1(3)) - (p3(2)-p1(2))*(p2(3)-p1(3)));
                ds2 = ((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)));
                cartesianPoints(noHeightIdx(i), 3) = p1(3) + dd1 / ds1 * (chp(2) - p1(2)) ...
                                                        - dd2 / ds2 * (chp(1) - p1(1));
                break % Above needs y = Ax + b 
            end
       end    
    end
    
end

function newPoints = checkHorizontal()
end

function newPoints = checkVertical()
end

% Calculate complex parts for triangle equations ahead to save number of
% computations later
function triangleEquations = calculateTriangEq(triangArr, coords)
    % maybe also c1,c2,c3 parts
    masterEquations = zeros(length(triangArr), 4)
    cTriangleVariables = zeros(length(triangArr), 6)
    for i = 1:length(triangArr)
        cTriangleVariables(j,1) = (p2(1) - p1(1));
        cTriangleVariables(j,2) = (p2(2) - p1(2));
        cTriangleVariables(j,3) = (p3(1) - p2(1));
        cTriangleVariables(j,4) = (p3(2) - p2(2));
        cTriangleVariables(j,5) = (p1(1) - p3(1));
        cTriangleVariables(j,6) = (p1(2) - p3(2));

        p1 = coords(triangArr(j,1), :); % x,y 1
        p2 = coords(triangArr(j,2), :); % x,y 2
        p3 = coords(triangArr(j,3), :); % x,y 3
        masterEquations(j,1) = ((p2(1)-p1(1))*(p3(3)-p1(3)) - (p3(1)-p1(1))*(p2(3)-p1(3))); % function for all triangles to avoid repetative calculation
        masterEquations(j,2) = ((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)));
        masterEquations(j,3) = ((p2(2)-p1(2))*(p3(3)-p1(3)) - (p3(2)-p1(2))*(p2(3)-p1(3)));
        masterEquations(j,4) = ((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)));
    end
end

function probDist = evaluateZeroSumGame(pointsUtils)
    % evaluate game
end