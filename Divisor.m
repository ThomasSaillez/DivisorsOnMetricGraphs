classdef Divisor < handle
    % This is a class made to handle divisors on metric graphs in Matlab. 
    % 
    % The positions of the chips of the divisor are stored has a pair of
    % vectors. The first one contains the indices of the edges (in the 
    % underlying metric graph) where the chips lie and the second one is a 
    % double width vector containing the distance from the beginning of the 
    % edge to the chips and the distance from the chips to the end of the 
    % edge. A third vector contains the degrees of those chips.
    %
    % When a position can be written in multiple ways, the methods will
    % write it with the smallest number of chips possible and with the
    % lowest edge indices possible. The chips are sorted by increasing edge
    % indices and then by increasing distance
    %
    % Note that most methods will also accept having input which are not 
    % written as above.

    properties (Access = protected)
        metricGraph
        edgeIndexVector
        distanceVector
        degreeVector
    end

    methods (Static, Access = private)
        function segmentVector = turnToSegments(edgeVector, segmentToEdge, segmentLengths)
            segmentVector = zeros(1, width(segmentToEdge));
            for i = 1:width(edgeVector)
                segmentVector = segmentVector + (edgeVector(i)*(segmentToEdge == i)) .* segmentLengths;
            end    
        end
    end

    methods (Access = protected)
        function reduceVectors(obj)
            % This function will shorten the three vectors of the chip state
            % in a unique reduced way

            % First step is to take chips that are on a edge to the incident edge of lowest index    
            [obj.edgeIndexVector, obj.distanceVector] = obj.metricGraph.reducePosition(obj.edgeIndexVector, obj.distanceVector);
            % In second, we sort the vector first by edge index then by lengths
            % on the edge
            [~, sorted]=sort(obj.distanceVector(1, :));
            obj.degreeVector = obj.degreeVector(sorted);
            obj.distanceVector = obj.distanceVector(:, sorted);
            obj.edgeIndexVector = obj.edgeIndexVector(sorted);
            [~, sorted] = sort(obj.edgeIndexVector);
            obj.degreeVector = obj.degreeVector(sorted);
            obj.distanceVector = obj.distanceVector(:, sorted);
            obj.edgeIndexVector = obj.edgeIndexVector(sorted);
            % We have a sorted vector, note that two chips on the same position
            % are consecutive
            % So in third, we make all the degree of chips at the same position
            % go two only one place

            fusionCheck(obj);
            % Last, we delete all the chips with degree 0
            startLength = length(obj.degreeVector);
            for i = 1:startLength
                if obj.degreeVector(startLength - i + 1) == 0
                    removeChipProtected(obj, startLength - i + 1);
                end    
            end

        end

        function fusionCheck(obj)
            for i = 1:length(obj.edgeIndexVector)-1
                if (obj.edgeIndexVector(i) == obj.edgeIndexVector(i+1)) && (obj.distanceVector(1, i) == obj.distanceVector(1, i+1))
                    obj.degreeVector(i+1) = obj.degreeVector(i+1) + obj.degreeVector(i);
                    obj.degreeVector(i) = 0;
                end    
            end
        end      

       function addChipPrivate(obj, edgeIndex, distanceFromTail, degree)
            if edgeIndex > obj.metricGraph.getNumEdges
                error('The edge index has to be lower than the number of edges of the graph. The given index is %d while the number of edges is %d.',edgeIndex, obj.metricGraph.getNumEdges)
            end    
            if edgeIndex < 1
                error('The edge index has to be strictly positive. The given index is %d.', edgeIndex)
            end
            if distanceFromTail < 0
                error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTail)
            end
            if distanceFromTail > obj.metricGraph.getLength(edgeIndex)
                error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTail, obj.metricGraph.getLength(edgeIndex))
            end
            obj.degreeVector = [obj.degreeVector degree];
            obj.distanceVector = [obj.distanceVector [distanceFromTail; obj.metricGraph.getLength(edgeIndex) - distanceFromTail]];
            obj.edgeIndexVector = [obj.edgeIndexVector edgeIndex];
       end

       function [distance, distanceFromEnd] = getNextChip(obj, edgeIndex, length, towardsHead)
           if nargin == 3
               towardsHead = true;
           end    
           if ~isa(towardsHead, 'logical')
               error('The 4th input as to be a logical, given input is %s', class(towardsHead))
           end
           chipsOnEdge = getInterior(obj, edgeIndex);
           if towardsHead
               if isempty(chipsOnEdge)
                   distanceFromEnd = 0;
                   distance = obj.metricGraph.getLength(edgeIndex);
               else    
                   positionOfPoint = find(length >= obj.distanceVector(1, chipsOnEdge), 1, "last");
                   if width(positionOfPoint) ~= 0
                        if positionOfPoint < width(chipsOnEdge)
                            distance = obj.distanceVector(1, chipsOnEdge(positionOfPoint + 1));
                            distanceFromEnd = obj.distanceVector(2, chipsOnEdge(positionOfPoint + 1));
                        else
                            distanceFromEnd = 0;
                            distance = obj.metricGraph.getLength(edgeIndex);
                        end
                   else
                        distance = obj.distanceVector(1, chipsOnEdge(1));
                        distanceFromEnd = obj.distanceVector(2, chipsOnEdge(1));
                    end
               end
           else
               if isempty(chipsOnEdge)
                   distance = 0;
                   distanceFromEnd = obj.metricGraph.getLength(edgeIndex);
               else 
                    positionOfPoint = find(length <= obj.distanceVector(1, chipsOnEdge),1,"first");
                    if width(positionOfPoint) ~= 0
                        if positionOfPoint > 1
                            distance = obj.distanceVector(1, chipsOnEdge(positionOfPoint - 1));
                            distanceFromEnd = obj.distanceVector(2, chipsOnEdge(positionOfPoint - 1));
                        else
                            distance = 0;
                            distanceFromEnd =  obj.metricGraph.getLength(edgeIndex);
                        end
                    else   
                        distance = obj.distanceVector(1, chipsOnEdge(width(chipsOnEdge)));
                        distanceFromEnd = obj.distanceVector(2, chipsOnEdge(width(chipsOnEdge)));
                    end
               end    
           end    
       end

       function outputSegment = exploreSegment(obj, pointEdgeIndex, pointDistance, towardsHead)
           [nextStop, ~] = obj.getNextChip(pointEdgeIndex, pointDistance, towardsHead);
           outputSegment = Segment(obj.metricGraph, pointEdgeIndex, pointDistance, nextStop);
       end

       function removeChipProtected(obj, chipIndex)
            if chipIndex > width(obj.degreeVector)
                error('You cannot have a chip index higher than number of chips. Given index is %d and the number of chips is %d.', chipIndex, width(obj.degreeVector))
            end 
            if chipIndex < 1
                error('You cannot have a chip index lower than 1. Given index is %d.', chipIndex)
            end    
            obj.degreeVector(chipIndex)=[];
            obj.edgeIndexVector(chipIndex)=[];
            obj.distanceVector(:, chipIndex)=[];
        end
   end

   
%---------------------------------------------------------------------------------------------    
    
    methods
        % Constructor
        function obj = Divisor(metricGraph)
            % This is the constructor of the class Divisor.
            % D = Divisor(A) creates a divisor defined on the graph A.
            %
            % See also MetricGraph.

            if ~isa(metricGraph,'MetricGraph')
                error('The input has to be a metric graph, the input class is %s.', class(metricGraph))
            end
            if ~metricGraph.isConnected
                error('The input metric graph has to be connected and it is not.')
            end    
            obj.metricGraph = metricGraph;
            obj.edgeIndexVector = zeros(1,0);
            obj.distanceVector = zeros(2,0);
            obj.degreeVector = zeros(1,0);
        end
        
        function addChip(obj, edgeIndex, distanceFromTail, degree)
            % myDivisor.addChip(edgeIndex, distanceFromTail, degree) will
            % increase the divisor myDivisor by degree at the position
            % (edgeIndex, distanceFromTail). Then it will try to simplify
            % the storage of the divisor.
            %
            % See also removeChip.
            
            if height(distanceFromTail) == 2
                distanceFromTail = distanceFromTail(1,:);
            end    
            if edgeIndex > obj.metricGraph.getNumEdges
                error('The edge index has to be lower than the number of edges of the graph. The given index is %d while the number of edges is %d.',edgeIndex, obj.metricGraph.getNumEdges)
            end    
            if edgeIndex < 1
                error('The edge index has to be strictly positive. The given index is %d.', edgeIndex)
            end
            if distanceFromTail < 0
                error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTail)
            end
            if distanceFromTail > obj.metricGraph.getLength(edgeIndex)
                error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTail, obj.metricGraph.getLength(edgeIndex))
            end
            addChipPrivate(obj, edgeIndex, distanceFromTail, degree);
            reduceVectors(obj);
        end
             
        function removeChip(obj, chipIndex)
            % myDivisor.removeChip(chipIndex) removes the chip of index
            % chipIndex of myDivisor. The order of the chips is explained
            % in the class.
            %
            % See also addChip and show.

            if chipIndex > width(obj.degreeVector)
                error('You cannot have a chip index higher than number of chips. Given index is %d and the number of chips is %d.', chipIndex, width(obj.degreeVector))
            end 
            if chipIndex < 1
                error('You cannot have a chip index lower than 1. Given index is %d.', chipIndex)
            end    
            obj.degreeVector(chipIndex)=[];
            obj.edgeIndexVector(chipIndex)=[];
            obj.distanceVector(:, chipIndex)=[];
        end
        
        function show(obj)
            % myDivisor.show will display a matrix containing on the first
            % line the edge indices of myDivisor, then the second and third
            % line will contain respectively the distance from the
            % start/the end of the edge to the chip and the fourth line
            % will contain the degrees of the corresponding chips.
            %
            % See also plot, tikzPlot.

            disp([obj.edgeIndexVector; obj.distanceVector; obj.degreeVector])
        end
        
        function chipVector = getInterior(obj, edgeIndex)
            % myDivisor.getInterior(edgeIndex) output the indices of the
            % chips that lies in the interior of the edge edgeIndex.
            %
            % See also degreeOfPoint, getInteriorDistance.

            if edgeIndex > obj.metricGraph.getNumEdges
                error('You cannot have an edge index larger than the number of edges on the graph. Given index is %d while the number of edges is %d.', edgeIndex, obj.metricGraph.getNumEdges)
            end
            if edgeIndex < 1 
                error('You cannot have an edge index strictly smaller than 1. Given index is %d.', edgeIndex)
            end
            chipVector = intersect(intersect(find(obj.edgeIndexVector == edgeIndex),find(obj.distanceVector(1,:))),find(obj.distanceVector(2,:)));
        end
        
        function chipVector = getInteriorDistance(obj, edgeIndex, needSecondLine)
            % myDivisor.getInteriorDistance(edgeIndex, needSecondLine)
            % outputs a vector containing the distances of the chips that
            % lay inside the interior of the edge edgeIndex. If
            % needSecondLine is false, then it will be a simple width
            % vector, if true or omitted it will be a double width vector.
            %
            % See also getInterior.

            if nargin == 2
                needSecondLine = true;
            end 
            if ~isa(needSecondLine, 'logical')
                error('Third input either needs to be logical or empty. Given input is a %s', class(needSecondLine))
            end
            if edgeIndex > obj.metricGraph.getNumEdges
                error('You cannot have an edge index larger than the number of edges on the graph. Given index is %d while the number of edges is %d.', edgeIndex, obj.metricGraph.getNumEdges)
            end
            if edgeIndex < 1 
                error('You cannot have an edge index strictly smaller than 1. Given index is %d.', edgeIndex)
            end
            if needSecondLine
                chipVector = obj.distanceVector(:, getInterior(obj, edgeIndex));
            else    
                chipVector = obj.distanceVector(1, getInterior(obj, edgeIndex));
            end
        end

        function cloned = clone(obj)
            % myDivisor.clone outputs a copy of the divisor myDivisor.

            if ~isa(obj,'Divisor')
                error('The input has to be a Divisor, the input class is %s.', class(obj))
            end
            cloned = Divisor(obj.metricGraph);
            cloned.degreeVector = obj.degreeVector;
            cloned.distanceVector = obj.distanceVector;
            cloned.edgeIndexVector = obj.edgeIndexVector;
        end

        function areSame = sameMetricGraph(obj1, obj2)
            % sameMetricGraph(Div1, Div2) outputs true if Div1 has the same
            % divisor as Div2 and false otherwise.

            if ~isa(obj1,'Divisor')
                error('The first input has to be a Divisor, the input class is %s.', class(obj))
            end
            if ~isa(obj2,'Divisor')
                error('The second input has to be a Divisor, the input class is %s.', class(obj))
            end
            if obj1.metricGraph == obj2.metricGraph
                areSame = true;
            else
                areSame = false;
            end    
        end

        function sum = addition(obj1, obj2)
            % addition(div1, div2) outputs a divisor which is the sum of
            % div1 and div2.
            %
            % See also scalarMultiplication.

            if ~isa(obj1,'Divisor')
                error('The first input has to be a Divisor, the input class is %s.', class(obj))
            end
            if ~isa(obj2,'Divisor')
                error('The second input has to be a Divisor, the input class is %s.', class(obj))
            end
            if ~sameMetricGraph(obj1, obj2)
                error('The input Divisors have to be defined on the same MetricGraph.')
            end    
            sum = obj1.clone;
            sum.addChip(obj2.edgeIndexVector,obj2.distanceVector(1,:),obj2.degreeVector);
        end

        function scalarMultiplication(obj, scalar)
            % myDivisor.scalarMultiplication(scalar) multiplies the degrees
            % of myDivisor by scalar.
            %
            % See also addition.

            if scalar ~= round(scalar)
                error('The second input has to be integer. Given input is %s.', scalar)
            end    
            if scalar == 0
                obj.degreeVector = [];
                obj.edgeIndexVector = [];
                obj.distanceVector = [];
            else
                obj.degreeVector = scalar * obj.degreeVector;
            end    
        end
        
        function degree = degreeOfDivisor(obj)
            % myDivisor.degreeOfDivisor outputs the degree of the divisor,
            % which is the total sum of the degrees of the individual chips
            % of the divisor.
            %
            % See also degreeOfPoint.

            degree = sum(obj.degreeVector);
        end

        function degree = degreeOfPoint(obj, edgeIndex, distance)
            % myDivisor.degreeOfPoint(edgeIndex, distance) outputs the
            % value of myDivisor at the point described by the pair
            % (edgeIndex, distance).
            %
            % See also degreeOfDivisor.
            
            if edgeIndex > obj.metricGraph.getNumEdges
                error('You cannot have an edge index larger than the number of edges on the graph. Given index is %d while the number of edges is %d.', edgeIndex, obj.metricGraph.getNumEdges)
            end
            if edgeIndex < 1 
                error('You cannot have an edge index strictly smaller than 1. Given index is %d.', edgeIndex)
            end
            if distance(1) < 0
                error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTail)
            end
            if distance(1) > obj.metricGraph.getLength(edgeIndex)
                error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTail, obj.metricGraph.getLength(edgeIndex))
            end

            degree = 0;
            if height(distance) == 1
                distance = [distance ; obj.metricGraph.getLength(edgeIndex)-distance];
            end    
            [edgeIndex, distance] = obj.metricGraph.reducePosition(edgeIndex, distance);
            matchingEdge = find(obj.edgeIndexVector == edgeIndex);
            matchingDistance = find(obj.distanceVector(1,:) == distance(1));
            matching = intersect(matchingDistance, matchingEdge);
            if ~isempty(matching)
                degree = obj.degreeVector(matching);
            end    
        end

        function toExplore = leavingSegments(obj, edgeIndex, distance)
            % myDivisor.leavingSegments(edgeIndex, distance) outputs all the
            % segments that start from the point described by the pair
            % (edgeIndex, distance).
            %
            % See also Segment.
             if edgeIndex > obj.metricGraph.getNumEdges
                error('You cannot have an edge index larger than the number of edges on the graph. Given index is %d while the number of edges is %d.', edgeIndex, obj.metricGraph.getNumEdges)
            end
            if edgeIndex < 1 
                error('You cannot have an edge index strictly smaller than 1. Given index is %d.', edgeIndex)
            end
            if distance(1) < 0
                error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTail)
            end
            if distance(1) > obj.metricGraph.getLength(edgeIndex)
                error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTail, obj.metricGraph.getLength(edgeIndex))
            end

            if height(distance) == 1
                distance = [distance ; obj.metricGraph.getLength(edgeIndex)-distance];
            end
            toExplore = zeros(0);
            if distance(1) == 0
                listOfIncidents = obj.metricGraph.getIncidents(obj.metricGraph.getTail(edgeIndex));
                for i = 1:width(listOfIncidents)
                    if obj.metricGraph.getIncidenceMatrix(obj.metricGraph.getTail(edgeIndex), listOfIncidents(i))==-1
                        D = exploreSegment(obj, listOfIncidents(i), 0, true);
                    else
                        D = exploreSegment(obj, listOfIncidents(i), obj.metricGraph.getLength(listOfIncidents(i)), false);
                    end    
                    toExplore = [toExplore clone(D)];%#ok
                end
            elseif distance(2) == 0
                listOfIncidents = obj.metricGraph.getIncidents(obj.metricGraph.getHead(edgeIndex));
                for i = 1:width(listOfIncidents)
                    if obj.metricGraph.getIncidenceMatrix(obj.metricGraph.getHead(edgeIndex),listOfIncidents(i))==-1
                        D = exploreSegment(obj, listOfIncidents(i), 0, true);
                    else
                        D = exploreSegment(obj, listOfIncidents(i), obj.metricGraph.getLength(listOfIncidents(i)), false);
                    end    
                    toExplore = [toExplore clone(D)];%#ok
                end
            else
                D = exploreSegment(obj, edgeIndex, distance(1), false);
                toExplore = [toExplore clone(D)];
                D = exploreSegment(obj, edgeIndex, distance(1), true);
                toExplore = [toExplore clone(D)];
           end
        end

        function dharEffective(obj, edgeIndex, distance, X, Y, outputVisual)
            % myDivisor.dharEffective(edgeIndex, distance) applies Dhar's
            % algorithm to myDivisor at the point described by the pair
            % (edgeIndex, distance).
            %
            % myDivisor.dharEffective(edgeIndex, distance, X, Y) does the
            % same and plots the steps of the algorithm.
            %
            % myDivisor.dharEffective(edgeIndex, distance, X, Y, outputVisual) 
            % does the same but if outputVisual is false then it uses
            % tikzPlot instead of plot.
            %
            % See also plot, tikzPlot.

            if nargin == 5
                outputVisual = true;
            end    
            if obj.metricGraph.isMultigraph
                error('Can not apply Dhar algorithm on a multigraph.')
            end    
            if sum(obj.degreeVector < 0) ~= 0
                error('The input divisor has to be effective. The given divisor is not.')
            end
             if edgeIndex > obj.metricGraph.getNumEdges
                error('You cannot have an edge index larger than the number of edges on the graph. Given index is %d while the number of edges is %d.', edgeIndex, obj.metricGraph.getNumEdges)
            end
            if edgeIndex < 1 
                error('You cannot have an edge index strictly smaller than 1. Given index is %d.', edgeIndex)
            end
            if distance(1) < 0
                error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTail)
            end
            if distance(1) > obj.metricGraph.getLength(edgeIndex)
                error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTail, obj.metricGraph.getLength(edgeIndex))
            end

            if height(distance) == 1
                distance = [distance; obj.metricGraph.getLength(edgeIndex) - distance];
            end
            [edgeIndex, distance] = reducePosition(obj.metricGraph, edgeIndex, distance);
            % We have two arrays for Dhar's algorithm. We call segment
            % a minimal subset of an edge bounded by points of the divisors
            % and/or vertices. The first array "toExplore" contains the oriented segments we are
            % about to explore. The second "toFire" contains the oriented segments that
            % stopped "Dhar's fire", and is used to know which segment should not be explored.
          
            obj.addChip(edgeIndex, distance(1), 1)
            while true
                doNotExplore = zeros(0);
                toFireAsDiv = Divisor(obj.metricGraph);
                toExplore = obj.leavingSegments(edgeIndex, distance);
                while width(toExplore) ~= 0
                    doNotExplore = [doNotExplore toExplore(1)]; %#ok
                    [edgeExplored, distanceExplored] = toExplore(1).getEnd;
                    if toFireAsDiv.degreeOfPoint(edgeExplored, distanceExplored) == -obj.degreeOfPoint(edgeExplored, distanceExplored) || obj.degreeOfPoint(edgeExplored,distanceExplored) == 0
                        if width(toExplore) >= 2
                            [endExploredEdge, endExploredDistance] = toExplore(1).getEnd;
                            for j = width(toExplore):-1:2
                                [endExploredEdgeLoop, endExploredDistanceLoop] = toExplore(j).getEnd;
                                if endExploredEdge == endExploredEdgeLoop && endExploredDistance(1) == endExploredDistanceLoop(1)
                                    doNotExplore = [doNotExplore toExplore(j)]; %#ok
                                    toExplore(j) = [];
                                end    
                            end
                        end
                        toBeChecked = leavingSegments(obj, edgeExplored, distanceExplored);
                        for k = 1:width(toBeChecked)
                            isNotExplored = true;
                            for j = 1:width(doNotExplore)
                                if areSameUnoriented(doNotExplore(j), toBeChecked(k), obj.metricGraph)
                                    isNotExplored = false;
                                    invert(doNotExplore(j))
                                    toUnfireAsDiv = turnToDivisor(doNotExplore(j), obj.metricGraph);
                                    toFireAsDiv = addition(toFireAsDiv, toUnfireAsDiv);
                                    delete(toUnfireAsDiv);
                                    doNotExplore(j) = []; %A segment can only be reached twice
                                    break
                                end   
                            end
                            if isNotExplored && width(toExplore) >= 2
                                for j = 2:width(toExplore)
                                    if areSameUnoriented(toExplore(j), toBeChecked(k), obj.metricGraph)
                                        isNotExplored = false;
                                        toExplore(j) = []; %A segment can only be reached twice
                                        break
                                    end   
                                end
                            end    
                            if isNotExplored
                               toExplore = [toExplore toBeChecked(k)];%#ok
                            end 
                        end   
                    else
                        invert(toExplore(1));
                        toExploreAsDiv = turnToDivisor(toExplore(1), obj.metricGraph);
                        result = addition(toFireAsDiv, toExploreAsDiv);
                        toFireAsDiv.degreeVector = result.degreeVector;
                        toFireAsDiv.distanceVector = result.distanceVector;
                        toFireAsDiv.edgeIndexVector = result.edgeIndexVector;
                        delete(result)
                        delete(toExploreAsDiv);
                    end
                    toExplore(1)=[];
                end
                if width(doNotExplore) == 0
                    break
                else    
                    delete(toFireAsDiv)
                    firingLength = zeros(1,width(doNotExplore));
                    for i = 1:width(doNotExplore)
                        firingLength(i) = lengthToStop(doNotExplore(i), edgeIndex, distance);
                    end
                    divisorToAdd = Divisor(obj.metricGraph);
                    for i = 1:width(doNotExplore)
                        changeLengthTo(doNotExplore(i), min(firingLength));
                        divisorToAdd = addition(divisorToAdd, turnToDivisor(doNotExplore(i),obj.metricGraph));
                    end
                    result = addition(obj, divisorToAdd);
                    obj.degreeVector = result.degreeVector;
                    obj.distanceVector = result.distanceVector;
                    obj.edgeIndexVector = result.edgeIndexVector;
                    delete(result)
                    if nargin >= 5
                        obj.addChip(edgeIndex, distance(1), -1)
                        if outputVisual
                            obj.plot(X, Y)
                            pause(1)
                        else
                            obj.tikzPlot(X, Y)
                        end    
                        obj.addChip(edgeIndex, distance(1), 1)
                    end
                end    
            end
        obj.addChip(edgeIndex,distance(1),-1)    
        end

        function plot(obj, X, Y)
            % myDivisor.plot(X, Y) creates a visual plot of myDivisor where 
            % the i-th vertex of the underlying metric graph is placed at 
            % coordinates (X(i), Y(i)).
            %
            % myDivisor.plot executes myDivisor.plot(X, Y) with X and Y
            % being the laplacian coordinates of the underlying metric 
            % graph.
            %
            % See also tikzPlot.

            if nargin == 1
                [X,Y] = obj.metricGraph.getLaplacianCoordinates;
            end
            if obj.metricGraph.getNumVertices ~= width(X)
                error('The width of the vector X has to be the number of vertices of the underlying graph. The graph has %d vertices and the vector has width %s.',obj.metricGraph.getNumVertices, width(X))
            end
            if obj.metricGraph.getNumVertices ~= width(Y)
                error('The width of the vector Y has to be the number of vertices of the underlying graph. The graph has %d vertices and the vector has width %s.',obj.metricGraph.getNumVertices, width(Y))
            end
            plot(obj.metricGraph, X, Y)
            hold on;
            for i = 1:width(obj.degreeVector)
                fromVertex = find(obj.metricGraph.getIncidenceMatrix(':', obj.edgeIndexVector(i)) == -1);
                toVertex = find(obj.metricGraph.getIncidenceMatrix(':', obj.edgeIndexVector(i)) == 1);
                proportion = obj.distanceVector(1,i)/sum(obj.distanceVector(:,i));
                XCoordinate = X(fromVertex) + proportion*(X(toVertex)-X(fromVertex));
                YCoordinate = Y(fromVertex) + proportion*(Y(toVertex)-Y(fromVertex));
                degreeString = int2str (obj.degreeVector(i));
                plot(XCoordinate, YCoordinate, '.', 'MarkerSize', 50, 'Color',[0.3 0.3 0.3]);
                text(XCoordinate, YCoordinate, degreeString,'color','white',HorizontalAlignment='center')
            end
            hold off;
        end

        function tikzPlot(obj, X, Y, showOne)
            % myDivisor.tikzPlot(X, Y, showOne) creates the Tikz 
            % code that produces the visual plot of the graph myGraph where
            % the i-th vertex is located at coordinates (X(i), Y(i)).
            %
            % myDivisor.tikzPlot executes myDivisor.tikzPlot(X, Y) with X 
            % and Y being the laplacian coordinates.
            % 
            % One can ommit the input showOne. If this input is set 
            % to false, the code doesn't print the degree of the points
            % when they are 1.
            %
            % The Tikz code is put in the file named "output.txt", in your
            % current folder. If the file doesn't exist then it is created.
            %
            % See also plot.

            if nargin == 1
                [X,Y] = obj.metricGraph.getLaplacianCoordinates;
            end
            if obj.metricGraph.getNumVertices ~= width(X)
                error('The width of the vector X has to be the number of vertices of the underlying graph. The graph has %d vertices and the vector has width %s.',obj.metricGraph.getNumVertices, width(X))
            end
            if obj.metricGraph.getNumVertices ~= width(Y)
                error('The width of the vector Y has to be the number of vertices of the underlying graph. The graph has %d vertices and the vector has width %s.',obj.metricGraph.getNumVertices, width(Y))
            end
            if nargin <= 3
                showOne = false;
            end
            if ~isa(showOne, 'logical')
               error('The fourth input has either to be empty or logical. The input class is %s.', class(needLastTikzLine))
            end   
            tikzPlot(obj.metricGraph, X, Y, false)
            fileID = fopen("output.txt",'a+');
            for i = 1:width(obj.degreeVector)
                fromVertex = find(obj.metricGraph.getIncidenceMatrix(':',obj.edgeIndexVector(i)) == -1);
                toVertex = find(obj.metricGraph.getIncidenceMatrix(':',obj.edgeIndexVector(i)) == 1);
                proportion = obj.distanceVector(1,i)/sum(obj.distanceVector(:,i));
                coordinateX = X(fromVertex) + proportion*(X(toVertex)-X(fromVertex));
                coordinateY = Y(fromVertex) + proportion*(Y(toVertex)-Y(fromVertex));
                degreeString = int2str (obj.degreeVector(i));
                if showOne || obj.degreeVector(i) ~= 1
                    fprintf(fileID, '    \\filldraw[black!70, thick] (%d, %d) circle (4pt) node[white, scale = 0.8]{%s}; \n', coordinateX, coordinateY, degreeString);
                else
                    fprintf(fileID, '    \\filldraw[black!70, thick] (%d, %d) circle (4pt); \n', coordinateX, coordinateY);
                end    
            end
            fprintf(fileID, '\\end{tikzpicture}\n');
            fclose(fileID);
        end
        
        function norm = pNorm(obj, p)
            % myDivisor.pNorm computes the Euclidean metric of the
            % corresponding divisor.
            %
            % myDivisor.pNorm(p) computes the same metric but instead of
            % taking the corresponding integral in L_2, it takes it in L_p.
            
            if obj.degreeOfDivisor ~= 0
                error('The divisor has to have degree 0. The given divisor has degree %i.', obj.degreeOfDivisor)
            end    
            if nargin == 1
                p = 2;
            end
            [slopeVector, segmentLengths] = obj.functionFromDivisor;
            norm = (sum(segmentLengths*(abs(slopeVector).^p))).^(1/p);
        end

        function areEquivalent = areEquivalent(obj1, obj2)
            % areEquivalent(div1, div2) tests if div1 and div2 are the same
            % divisors when seen in the Picard group.

            tolerance = 10^(-15);
            if ~isa(obj1,'Divisor')
                error('The first input has to be a Divisor, the input class is %s.', class(obj))
            end
            if ~isa(obj2,'Divisor')
                error('The second input has to be a Divisor, the input class is %s.', class(obj))
            end
            if ~sameMetricGraph(obj1, obj2)
                error('The input Divisors have to be defined on the same MetricGraph.')
            end
            if obj1.degreeOfDivisor == obj2.degreeOfDivisor
                difference = obj2.clone;
                difference.scalarMultiplication(-1);
                difference = addition(difference, obj1);
                [slopeVector, ~] = difference.functionFromDivisor;
                areEquivalent = all((slopeVector - round(slopeVector)) < tolerance);
            else
                areEquivalent = false;
            end    
        end

        function [slopeVector, segmentLengths] = functionFromDivisor(obj)
            % myDivisor.functionFromDivisor outputs the vector containing
            % all the slopes that a function would have if it was a
            % principal divisor. Slopes are oriented the same way as the
            % edges they lie in and are ordered from first edge to last
            % edge from the start of an edge to its end.
            

            if obj.degreeOfDivisor ~= 0
                error('The input divisor has degree nonzero, given degree is %i.', obj.degreeOfDivisor)
            end    
            supportVector = ones(1,obj.metricGraph.getNumEdges);
            segmentLengths = [];
            for i = 1:obj.metricGraph.getNumEdges
                interiorLengths = obj.getInteriorDistance(i, false);
                switch width(interiorLengths)
                    case 0
                        segmentLengths = [segmentLengths obj.metricGraph.getLength(i)]; %#ok
                    case 1
                        segmentLengths = [segmentLengths interiorLengths(1) obj.metricGraph.getLength(i) - interiorLengths(end)]; %#ok
                    otherwise
                        segmentLengths = [segmentLengths interiorLengths(1) diff(interiorLengths) obj.metricGraph.getLength(i) - interiorLengths(end)]; %#ok
                end    
                supportVector(i) = width(interiorLengths) + 1;
            end    
            matrixDimension = sum(supportVector);
            segmentToEdge = zeros(1, matrixDimension);
            pointerBis = 1;
            for i = 1:obj.metricGraph.getNumEdges
                segmentToEdge(pointerBis:pointerBis + supportVector(i)-1) = i;
                pointerBis = pointerBis + supportVector(i);
            end
            pointer = 1;
            slopeSystemMatrix = zeros(matrixDimension);
            slopeSystemVector = zeros(matrixDimension,1);
            pathMatrix = zeros(obj.metricGraph.getNumVertices, obj.metricGraph.getNumEdges);
            toExplore = obj.metricGraph.getIncidents(1);
            toExploreOrigins = ones(1,width(toExplore));
            explored = toExplore;
            while width(toExplore) ~= 0
                currentVertex = obj.metricGraph.getOtherVertex(toExplore(1), toExploreOrigins(1));
                listOfAdjacents = obj.metricGraph.getIncidents(currentVertex);
                listOfAdjacents = setdiff(listOfAdjacents, intersect(listOfAdjacents, explored));
                toExplore = [toExplore listOfAdjacents]; %#ok
                explored = [explored listOfAdjacents]; %#ok
                toExploreOrigins = [toExploreOrigins currentVertex*ones(1,width(listOfAdjacents))]; %#ok               
                % The next if will add the information about the vertex to
                % our matrices containing the paths.
                if any(pathMatrix(currentVertex,:))
                    loopExplored = pathMatrix(toExploreOrigins(1),:) - pathMatrix(currentVertex,:);
                    if obj.metricGraph.getTail(toExplore(1)) == toExploreOrigins(1)
                        loopExplored(toExplore(1)) = 1;
                    else
                        loopExplored(toExplore(1)) = -1;
                    end
                    slopeSystemMatrix(pointer, :) = Divisor.turnToSegments(loopExplored, segmentToEdge, segmentLengths);
                    pointer = pointer + 1;
                else
                    pathMatrix(currentVertex,:) = pathMatrix(toExploreOrigins(1),:);
                    if obj.metricGraph.getTail(toExplore(1)) == toExploreOrigins(1)
                        pathMatrix(currentVertex, toExplore(1)) = 1;
                    else
                        pathMatrix(currentVertex, toExplore(1)) = -1;
                    end    
                end
                toExplore(1) = [];
                toExploreOrigins(1) = [];
            end
            for i = 1:obj.metricGraph.getNumVertices - 1
                [edge, distance] = obj.metricGraph.vertexToPoint(i);
                slopeSystemVector(pointer) = obj.degreeOfPoint(edge, distance);
                headIncidents = obj.metricGraph.getHeadIncidents(i);
                for j = 1:width(headIncidents)
                    slopeSystemMatrix(pointer, find(segmentToEdge == headIncidents(j),1,"last")) = 1;
                end
                tailIncidents = obj.metricGraph.getTailIncidents(i);
                for j = 1:width(tailIncidents)
                    slopeSystemMatrix(pointer, find(segmentToEdge == tailIncidents(j),1,"first")) = -1;
                end 
                pointer = pointer + 1;
            end
            for i = 1:obj.metricGraph.getNumEdges
                interiorLengths = obj.getInteriorDistance(i, false);
                startPosition = find(segmentToEdge == i, 1) - 1;
                for j = 1:width(interiorLengths)
                    slopeSystemMatrix(pointer, startPosition + j) = 1; 
                    slopeSystemMatrix(pointer, startPosition + j + 1) = -1;
                    slopeSystemVector(pointer) = obj.degreeOfPoint(i,interiorLengths(j));
                    pointer = pointer + 1;
                end    
            end
            slopeVector = linsolve (slopeSystemMatrix, slopeSystemVector);
        end
    end
end