close all
format short
clear

runThesis();

% Function encapsulating outputs to keep Workspace empty
function runThesis()
    str = input("Do you want to change default setting? Y or enter: ",'s');
    if ~isempty(str)
        newN = input("Set number of points: "); % check if integer
        newSmooth = input("Set smoothness of points: ");
        newRuns = input("Set number of runs: "); % don't have to set all
        newHeight = input("Set possible max height of points: ");
        [n, smoothVal, numRuns, height] = setVariables(newN, newSmooth, newRuns, newHeight);
        clear newN newSmooth newRuns; % was in local scope
    else
        [n, smoothVal, numRuns, height] = setVariables(10, 0.2, 5, 500);
    end
    visuals = input("Do you want to see all outputs? Y or enter: ",'s');
    if ~isempty(visuals)
        visualsIdx = 1;
    else
        visualsIdx = numRuns;
    end
    
    [pointsArr, heightList] = generateEnvironment(smoothVal, n, height);
    [r, c] = find(pointsArr > -1);
    bools = pointsArr > -1;
    coords = [r, c, pointsArr(bools)]; % x,y,z
    
    triangArr = createTriangulations(r, c);
    [lineList, lLwCoords] = createListOfLines(coords, triangArr);
    sLWCoords = sortrows(lLwCoords, 'ascend');
    coordsForEqs = coords;
    [cTriangleVariables, triangleEquations] = calculateTriangEq(triangArr, coordsForEqs);

    for i = 1:numRuns
        if i ~= 1
            coords = cartesianPoints;
        end
        tic %TIMER START
        smootherpoints = smoothenPlane(coords, sLWCoords);
        cartesianPoints = makeCartesianProduct(smootherpoints);
        cartesianPoints = [cartesianPoints, -ones(size(cartesianPoints,1), 1)]; % x y -1
        [~, indx] = intersect(cartesianPoints(:,1:2), smootherpoints(:,1:2), 'rows'); % assign height to points already counted
        cartesianPoints(indx, 3) = smootherpoints(:, 3);
        toc %TIMER STOPS
    
        pointsPrint = ['Number of points in domain is: ', num2str(size(cartesianPoints,1))];
        disp(pointsPrint);
        currRun = ['Currently we are at iteration number ', num2str(i), ' out of ', num2str(numRuns), '.'];
        disp(currRun);
    
        cartesianPoints = calculateCPHeights(triangArr, coordsForEqs, cartesianPoints, triangleEquations, cTriangleVariables);
        [probDist, probDistt] = evaluateZeroSumGame(cartesianPoints);
        if visualsIdx == i
            displayTriangulated(triangArr, heightList, r, c, cartesianPoints, probDist, probDistt);
            visualsIdx = visualsIdx + 1;
            if i ~= numRuns
                disp("Press a key to continue!");
                pause;
            end
        end    
    end
    disp("Finished");
end

% Sets required variables and checks valid values
function [n, smoothVal, numRuns, height]  = setVariables(newN, newSmooth, newRuns, newHeight)
    n = checkValidity(newN, 4);
    smoothVal = checkValidity(newSmooth, 0.1);
    numRuns = checkValidity(newRuns, 1);
    height = checkValidity(newHeight, 0);
end

% Check input validity
function checkedValue = checkValidity(valueWanted, threshold)
    if valueWanted < threshold
       tooLow = ["Value too low setting to ", threshold];
       disp(tooLow);
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
    idx = sort([noCornerIdx(randperm(numel(noCornerIdx), n - 4)), cornerIdx]);
    heightRange = height + 1;
    heightList = randperm(heightRange, n) - 1;
    pointsArr(idx) = heightList;
end

% Returns delaunay triangulation using delaunay library
function triangArr = createTriangulations(r, c)
    triangArr = delaunay(r, c);
end

% Displays either 2D, 3D, or both
function displayTriangulated(triangArr, heightList, r, c, cartesianPoints, probDist, probDistt)
    str = input("Do you want 2D, 3D, or both visualizations? 2,3, or enter: ", 's');
    if (strcmp(str,"2"))
        figure('Name','2D Visualization','NumberTitle','off')
        triplot(triangArr, r, c);
        hold on
        plot(cartesianPoints(:, 1), cartesianPoints(:, 2), 'r.')
        xlabel("x coordinate")
        ylabel("y coordinate")
        hold off
    elseif (strcmp(str,"3"))    
        figure('Name','3D Visualization','NumberTitle','off')
        trimesh(triangArr, r, c, heightList);
        hold on
        plot3(cartesianPoints(:, 1), cartesianPoints(:, 2), cartesianPoints(:, 3), 'r.')
        xlabel("x coordinate")
        ylabel("y coordinate")
        zlabel("z coordinate")
        hold off
    else
        figure('Name','2D Visualization','NumberTitle','off')
        triplot(triangArr, r, c);
        hold on
        plot(cartesianPoints(:, 1), cartesianPoints(:, 2), 'r.')
        xlabel("x coordinate")
        ylabel("y coordinate")
        hold off
        figure('Name','3D Visualization','NumberTitle','off')
        trimesh(triangArr, r, c, heightList);
        hold on
        plot3(cartesianPoints(:, 1), cartesianPoints(:, 2), cartesianPoints(:, 3), 'r.')
        xlabel("x coordinate")
        ylabel("y coordinate")
        zlabel("z coordinate")
        hold off
    figure('Name','Player Probability Distribution','NumberTitle','off')
    stairs(probDist, 'LineWidth', 2, 'Marker', 'd', 'MarkerFaceColor', 'c');
    hold on
    stairs(probDistt, 'LineWidth', 2, 'Marker', 'd', 'MarkerFaceColor', 'r');
    hold off
    end
end

% Function returns lineList = list of lines with coords indices and
% lLwCoords = array of line coordinations (x1, y1, x2, y2)
function [lineList, lLwCoords] = createListOfLines(coords, triangArr)
    nTriangRows = size(triangArr, 1);
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
    lLwCoords = zeros(nOLines(1), 6);
    for i = 1:nOLines
        lLwCoords(i, :) = [coords(lineList(i,1), :), coords(lineList(i,2), :)]; %, coords(lineList(i,3), :)];
    end
end

% Find new intersections in both directions
function smootherpoints = smoothenPlane(coords, sLWCoords)
    smootherpoints = coords;
    for i = 1:size(coords, 1)
        verticalInt = sLWCoords((sLWCoords(:,1)<coords(i,1) & sLWCoords(:,4)>coords(i,1)) | ...
                        (sLWCoords(:,4)<coords(i,1) & sLWCoords(:,1)>coords(i,1)), :);

        for j = 1:size(verticalInt, 1)
            if i>1 && coords(i-1,1) == coords(i,1)
                break
            end
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
            if i>1 && coords(i-1,2) == coords(i,2)
                break
            end
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
    [~,ia,~] = unique(smootherpoints(:,1:2),'rows');
    smootherpoints = smootherpoints(ia, :);
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

% Using pre-calculated equations find triangle belonging to each unassigned
% point and assign them height
function cartesianPoints = calculateCPHeights(triangArr, coords, cartesianPoints, triangleEquations, cTriangleVariables)
    noHeightIdx = find(cartesianPoints(:,3) == -1);
    for i = 1:length(noHeightIdx)
       for j = 1:length(triangArr)
            p1 = coords(triangArr(j,1), :); % x,y 1
            p2 = coords(triangArr(j,2), :); % x,y 2
            p3 = coords(triangArr(j,3), :); % x,y 3
            chp = cartesianPoints(noHeightIdx(i), 1:2); % x,y checked
            c1 = cTriangleVariables(j,1) * (chp(2) - p1(2)) - cTriangleVariables(j,2) * (chp(1) - p1(1));
            c2 = cTriangleVariables(j,3) * (chp(2) - p2(2)) - cTriangleVariables(j,4) * (chp(1) - p2(1));
            c3 = cTriangleVariables(j,5) * (chp(2) - p3(2)) - cTriangleVariables(j,6) * (chp(1) - p3(1));
            if (c1 <= 0 && c2 <= 0 && c3 <= 0) || (c1 >= 0 && c2 >= 0 && c3 >= 0)
                cartesianPoints(noHeightIdx(i), 3) = p1(3) + triangleEquations(j,1) / triangleEquations(j,2) ...
                        * (chp(2) - p1(2)) - triangleEquations(j,3) / triangleEquations(j,4) * (chp(1) - p1(1));
                break
            end
       end    
    end
    
end

% Calculate complex parts for triangle equations ahead to save number of
% computations later
function [cTriangleVariables, triangleEquations] = calculateTriangEq(triangArr, coords)
    triangleEquations = zeros(length(triangArr), 4);
    cTriangleVariables = zeros(length(triangArr), 6);
    for i = 1:length(triangArr)
        p1 = coords(triangArr(i,1), :); % x,y 1
        p2 = coords(triangArr(i,2), :); % x,y 2
        p3 = coords(triangArr(i,3), :); % x,y 3        
        cTriangleVariables(i,1) = (p2(1) - p1(1));
        cTriangleVariables(i,2) = (p2(2) - p1(2));
        cTriangleVariables(i,3) = (p3(1) - p2(1));
        cTriangleVariables(i,4) = (p3(2) - p2(2));
        cTriangleVariables(i,5) = (p1(1) - p3(1));
        cTriangleVariables(i,6) = (p1(2) - p3(2));
        triangleEquations(i,1) = ((p2(1)-p1(1))*(p3(3)-p1(3)) - (p3(1)-p1(1))*(p2(3)-p1(3))); % function for all triangles to avoid repetative calculation
        triangleEquations(i,2) = ((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)));
        triangleEquations(i,3) = ((p2(2)-p1(2))*(p3(3)-p1(3)) - (p3(2)-p1(2))*(p2(3)-p1(3)));
        triangleEquations(i,4) = ((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)));
    end
end

% Calculate game probabilities
function [probDist, probDistt] = evaluateZeroSumGame(cartesianPoints)
    matrixSize = groupcounts(cartesianPoints(:, 2));
    occsX = matrixSize(1);
    occsY = size(cartesianPoints, 1) / occsX;
    A = reshape(cartesianPoints(:, 3), [occsY, occsX]);
    
    % Alice
    f = cat(1, zeros(occsX, 1), 1);
    Am = cat(2, A, -ones(occsY, 1));
    b = zeros(occsY, 1);
    Aeq = cat(1, ones(occsX, 1), 0)';
    beq = 1;
    lb = zeros(occsX, 1);
    ub = [];

    res = linprog(f, Am, b, Aeq, beq, lb, ub);
    probDist = res(1:occsX, 1);
    utilVal = res(end, 1);

    % Bob
    ft = cat(1, zeros(occsY, 1), 1); % 't' for two
    An = cat(2, transpose(A), -ones(occsX, 1));
    bt = zeros(occsX, 1);
    Aeqt = cat(1, ones(occsY, 1), 0)';
    beqt = 1;
    lbt = zeros(occsY, 1);
    ubt = [];

    rest = linprog(-ft, -An, bt, Aeqt, beqt, lbt, ubt);
    probDistt = rest(1:occsY, 1);
    utilValt = rest(end, 1);
end
