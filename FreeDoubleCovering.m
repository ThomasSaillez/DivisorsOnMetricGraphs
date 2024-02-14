classdef FreeDoubleCovering < MetricGraph & handle
    % This is a class made to handle free double coverings of metric graphs.
    % This class is inherited from MetricGraph, but one cannot modify the
    % graph with the same techniques.
    %
    % When constructing a double cover, the constructor will make two
    % copies of the base graph. The modifications you can do is make a pair
    % of edge cross between the two copies of the graph or stay in the same
    % copy.
    %
    % See also MetricGraph.

    properties (Access = private)
        baseGraph
    end
    properties (Access = protected)
        isCrossing
    end
    methods (Access = protected)
        function isConnected = isPreimageOfCycleConnected(obj, cycle)
            numberOfCrossings = cycle*obj.isCrossing;
            isConnected = mod(numberOfCrossings, 2);
        end

        function [firstComponent, liftedPathMatrix, parityVertices] = getLiftedTree(obj)
            [~, tree, pathMatrix] = obj.baseGraph.getBasisOfCycles;
            % relativeDistanceFromRoot can be negative, as going against
            % orientation decreases the distance, but it is a feature which
            % makes sure the correct edge is chosen when having both
            % orientations.
            relativeDistanceFromRoot = sum(transpose(pathMatrix));
            parityVertices = transpose(mod(pathMatrix * obj.isCrossing, 2));
            parityEdges = zeros(1, obj.getHalfNumEdges);
            for i = 1:obj.getHalfNumEdges
                vertices = obj.baseGraph.getIncidentVertices(i);
                [~, minIndex] = min(relativeDistanceFromRoot(vertices));
                parityEdges(i) = parityVertices(vertices(minIndex));
            end    
            firstComponent = [tree.*(1 - parityEdges) tree.*parityEdges];
            liftedPathMatrix = kron(eye(2), pathMatrix);
            for i = 1:obj.getHalfNumEdges
                if parityEdges(i)
                    temp = liftedPathMatrix(:, i);
                    liftedPathMatrix(:, i) = liftedPathMatrix(:, i + obj.getHalfNumEdges);
                    liftedPathMatrix(:, i + obj.getHalfNumEdges) = temp;
                end    
            end
            for i = 1:obj.getHalfNumVertices
                if parityVertices(i)
                    temp = liftedPathMatrix(i, :);
                    liftedPathMatrix(i, :) = liftedPathMatrix(i + obj.getHalfNumVertices, :);
                    liftedPathMatrix(i + obj.getHalfNumVertices, :) = temp;
                end    
            end
            parityVertices = [parityVertices 1-parityVertices];
        end    
    end

    methods
        %Constructor
        function obj = FreeDoubleCovering(metricGraph)
            % This is the constructor of FreeDoubleCovering.
            % myCover = FreeDoubleCovering(myGraph) makes a free double covering out
            % of the graph myGraph.
            %
            % See also MetricGraph.

            if ~isa(metricGraph,'MetricGraph')
                error('The input has to be a metric graph, the input class is %s.', class(metricGraph))
            end
            numVertices = 2*metricGraph.getNumVertices;
            obj@MetricGraph(numVertices);
            obj.baseGraph = metricGraph;            
            obj.incidenceMatrix = [metricGraph.getIncidenceMatrix zeros(metricGraph.getNumVertices, metricGraph.getNumEdges);zeros(metricGraph.getNumVertices, metricGraph.getNumEdges) metricGraph.getIncidenceMatrix];
            obj.lengths = [metricGraph.getLength metricGraph.getLength];
            obj.isCrossing = false(metricGraph.getNumEdges, 1);
            obj.isKnownConnected = true;
            obj.computedConnected = false;
            obj.isKnownMultigraph = metricGraph.isKnownMultigraph;
            obj.computedMultigraph = metricGraph.computedMultigraph;
        end
        
        function numVertices = getHalfNumVertices(obj)
            % myCover.getHalfNumVertices outputs half of the number of
            % vertices of myCover, which is the number of vertices of the
            % base graph.
            %
            % See also getNumEdges, getHalNumEdges.

            numVertices = obj.baseGraph.getNumVertices;
        end    
        
        function numEdges = getHalfNumEdges(obj)
            % myCover.getHalfNumEdges outputs half of the number of
            % edges of myCover, which is the number of edges of the base 
            % graph.
            %
            % See also getNumVertices, getHalNumVertices.

            numEdges = obj.baseGraph.getNumEdges;
        end    

        function makeEdgeCross(obj, edgeNumber)
            % myCover.makeEdgeCross(edgeNumber) makes that the edge
            % edgeNumber and the other copy of this edge go from one copy
            % of the base graph to the other.
            %
            % The first copy of a crossing edge will always go from the
            % first fiber to the second and the second copy will always go
            % from the second fiber to the first.
            %
            % See also makeEdgeUncross.

            if edgeNumber > obj.getNumEdges
                error('The edge number cannot be larger than the number of edges. The given number is %d.', edgeNumber)
            end  
            if edgeNumber < 1
               error('The edge number cannot be smaller than 1. The given number is %d.', edgeNumber)
            end
            obj.isKnownConnected = false;
            obj.isKnownMultigraph = false;
            if edgeNumber > obj.getHalfNumEdges
                edgeNumber = edgeNumber - obj.getHalfNumEdges;
            end
            if ~obj.isCrossing(edgeNumber)
                obj.isCrossing(edgeNumber) = true;
                vertex = find(obj.incidenceMatrix(:, edgeNumber) == 1);
                obj.incidenceMatrix(vertex, edgeNumber) = 0;
                obj.incidenceMatrix(vertex, edgeNumber + obj.getHalfNumEdges) = 1;
                obj.incidenceMatrix(vertex + obj.getHalfNumVertices, edgeNumber) = 1;
                obj.incidenceMatrix(vertex + obj.getHalfNumVertices, edgeNumber + obj.getHalfNumEdges) = 0;
            end    
        end

        function makeEdgeUncross(obj, edgeNumber)
            % myCover.makeEdgeUncross(edgeNumber) makes that the edge
            % edgeNumber and the other copy of this edge stay in the same
            % copy of the base graph.
            %
            % See also makeEdgeCross.

            if edgeNumber > obj.getNumEdges
                error('The edge number cannot be larger than the number of edges. The given number is %d.', edgeNumber)
            end  
            if edgeNumber < 1
               error('The edge number cannot be smaller than 1. The given number is %d.', edgeNumber)
            end
            obj.isKnownMultigraph = false;
            obj.isKnownConnected = false;
            if obj.isCrossing(edgeNumber)
                obj.isCrossing(edgeNumber) = false;
                vertex = find(obj.incidenceMatrix(:, edgeNumber + obj.getHalfNumEdges) == 1);
                obj.incidenceMatrix(vertex, edgeNumber) = 1;
                obj.incidenceMatrix(vertex, edgeNumber + obj.getHalfNumEdges) = 0;
                obj.incidenceMatrix(vertex + obj.getHalfNumVertices, edgeNumber) = 0;
                obj.incidenceMatrix(vertex + obj.getHalfNumVertices, edgeNumber + obj.getHalfNumEdges) = 1;
            end    
        end

        function [X, Y] = getLaplacianCoordinates(obj, scaleFactor)
            % [X, Y] = myCover.getLaplacianCoordinates(scaleFactor) outputs
            % coordiantes based on the second and third eigenvectors of the
            % laplacian matrix of the base graph of myCover.
            %
            % If scaleFactor is not precised, the default value used is 2.
            %
            % See also getLaplacianMatrix, plot, tikzPlot.

            if nargin < 2
                scaleFactor = 2;
            end
            [x, y] = getLaplacianCoordinates(obj.baseGraph, scaleFactor);
            y = y + 1.1*scaleFactor;
            X = [x -x];
            Y = [y -y];
        end

        function metricGraph = forgetCovering(obj)
            % myCover.forgetCovering outputs myCover has a MetricGraph.

            metricGraph = MetricGraph(obj.getNumVertices);
            metricGraph.incidenceMatrix = obj.incidenceMatrix;
            metricGraph.lengths = obj.lengths;
            metricGraph.isKnownConnected = obj.isKnownConnected;
            metricGraph.computedConnected = obj.computedConnected;
            metricGraph.isKnownMultigraph = obj.isKnownMultigraph;
            metricGraph.computedMultigraph = obj.computedMultigraph;
        end

        function [edgeIndex, distance] = involution(obj, edgeIndex, distance)
            % myCover.involution(edgeIndex, distance) computes the image of
            % the point described by the pair (edgeIndex, distance) under 
            % the special involution of myCover. If the distance is omitted
            % it is assummed to be 0.
            
            if edgeIndex > obj.getNumEdges
                error('The edge index has to be lower than the number of edges of the graph. The given index is %d while the number of edges is %d.', edgeIndex, obj.getNumEdges)
            end    
            if edgeIndex < 1
                error('The edge index has to be strictly positive. The given index is %d.', edgeIndex)
            end
            if nargin == 2
                distance = 0;
            end    
            if distance < 0
                error('The distance has to be positive. Given distance is %d.', distance)
            end
            if distance > obj.lengths(edgeIndex)
                error('The distance has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distance, obj.lengths(edgeIndex))
            end
            if edgeIndex > obj.getHalfNumEdges
                edgeIndex = edgeIndex - obj.getHalfNumEdges;
            else
                edgeIndex = edgeIndex + obj.getHalfNumEdges;
            end    
        end
        
        function vertex = involutionVertex(obj, vertex)
            % myCover.involutionVertex(vertexIndex) computes the image of
            % the given vertex by the special involution of myCover. 

            if vertex > obj.getNumVertices
                error('The vertex index has to be lower than the number of vertices of the graph. The given index is %d while the number of vertices is %d.', vertex, obj.getNumVertices)
            end    
            if vertex < 1
                error('The vertex index has to be strictly positive. The given index is %d.', vertex)
            end

            if vertex > obj.getHalfNumVertices
                vertex = vertex - obj.getHalfNumVertices;
            else
                vertex = vertex + obj.getHalfNumVertices;
            end
        end    

        function addEdge(obj) %#ok
            % One may not add an edge to a free double covering.
            %
            % See also forgetCovering.
            error('One may not add an edge to a free double covering.')
        end

        function removeEdge(obj) %#ok
            % One may not remove an edge from a free double covering.
            %
            % See also forgetCovering.
            error('One may not remove an edge from a free double covering.')
        end

        function removeAllEdges(obj) %#ok
            % One may not remove an edge from a free double covering.
            %
            % See also forgetCovering.
            error('One may not remove an edge from a free double covering.')
        end

        function removeAllEdgesSym(obj) %#ok
            % One may not remove an edge from a free double covering.
            %
            % See also forgetCovering.
            error('One may not remove an edge from a free double covering.')
        end
        
        function changeOrientation(obj) %#ok
            % One may not change an edge from a free double covering.
            %
            % See also forgetCovering.
            error('One may not change an edge from a free double covering.')
        end

        function cloned = clone(obj)
           % myCover.clone outputs a copy of the graph myGraph, which will
           % have exactly the same properties.

           if ~isa(obj,'FreeDoubleCovering')
                error('The input has to be a free double covering, the input class is %s.', class(obj))
           end    
           cloned = FreeDoubleCovering(obj.baseGraph);
           cloned.incidenceMatrix = obj.incidenceMatrix;
           cloned.lengths = obj.lengths;
           cloned.isCrossing = obj.isCrossing;
           cloned.isKnownConnected = obj.isKnownConnected;
           cloned.computedConnected = obj.computedConnected;
           cloned.isKnownLaplacianCoordinates = obj.isKnownLaplacianCoordinates;
           cloned.computedLaplacianCoordinatesX = obj.computedLaplacianCoordinatesX;
           cloned.computedLaplacianCoordinatesY = obj.computedLaplacianCoordinatesY;
        end

        function [prymMatrix, cycleMatrix] = getPrymVariety(obj)
            % [prymMatrix, cycleMatrix] = myCovering.getPrymVariety outputs
            % two informations.
            %
            % The prymMatrix is the Prym variety of
            % the covering.
            % 
            % The cycleMatrix outputs  a basis of the indenpendent cycles
            % of the kernel of the map sending cycles from the fibre to the
            % base graph. Every row will correspond to a cycle, every
            % column to an edge and a +1 tells that the edge is present in
            % the cycle in its orientation and a -1 that the edge is
            % present but with its opposite orientation.

            [firstComponent, liftedPathMatrix, firstComponentVertices] = obj.getLiftedTree;
            secondComponent = [firstComponent(obj.getHalfNumEdges + 1:obj.getNumEdges) firstComponent(1:obj.getHalfNumEdges)];
            bothComponentsEdges = firstComponent + secondComponent;
            edgesList = find(bothComponentsEdges == 0);
            cycleMatrix = zeros(length(edgesList)/2 - 1, obj.getNumEdges);
            firstCrossingEdge = [];
            i = 0;
            for j = 1:length(edgesList)/2
                i = i + 1 ;
                incidentVertices = obj.getIncidentVertices(edgesList(j));
                if firstComponentVertices(incidentVertices(1)) == firstComponentVertices(incidentVertices(2))
                    cycleMatrix(i, :) = liftedPathMatrix(incidentVertices(1), :) ...
                        - liftedPathMatrix(incidentVertices(2), :) ...
                        - liftedPathMatrix(obj.involutionVertex(incidentVertices(1)), :) ...
                        + liftedPathMatrix(obj.involutionVertex(incidentVertices(2)), :) ;
                    if incidentVertices(1) == obj.getTail(edgesList(j))
                        cycleMatrix(i, edgesList(j)) = 1;
                        cycleMatrix(i, edgesList(j) + obj.getHalfNumEdges) = -1;
                    else
                        cycleMatrix(i, edgesList(j)) = -1;
                        cycleMatrix(i, edgesList(j) + obj.getHalfNumEdges) = 1;
                    end
                else
                    if isempty(firstCrossingEdge)
                        firstCrossingEdge = edgesList(j);
                        firstCrossingVertices = [obj.getTail(firstCrossingEdge) obj.getHead(firstCrossingEdge)];
                        if firstComponentVertices(firstCrossingVertices(1))
                            firstCrossingOrientation = 1;
                        else
                            firstCrossingOrientation = -1;
                            temp = firstCrossingVertices(1);
                            firstCrossingVertices(1) = firstCrossingVertices(2);
                            firstCrossingVertices(2) = temp;
                        end
                        i = i - 1;
                    else
                        cycleMatrix(i, :) = + liftedPathMatrix(obj.involutionVertex(firstCrossingVertices(2)), :) ...
                                - liftedPathMatrix(obj.involutionVertex(firstCrossingVertices(1)), :) ...
                                - liftedPathMatrix(firstCrossingVertices(2), :) ...
                                + liftedPathMatrix(firstCrossingVertices(1), :);
                        if firstComponentVertices(obj.getTail(edgesList(j)))
                            cycleMatrix(i, :) = cycleMatrix(i, :) ...
                                - liftedPathMatrix(incidentVertices(1), :) ...                                
                                - liftedPathMatrix(obj.involutionVertex(incidentVertices(2)), :) ...
                                + liftedPathMatrix(obj.involutionVertex(incidentVertices(1)), :) ...
                                + liftedPathMatrix(incidentVertices(2), :) ;
                            cycleMatrix(i, edgesList(j)) = 1;
                            cycleMatrix(i, edgesList(j) + obj.getHalfNumEdges) = -1;
                        else
                             cycleMatrix(i, :) = cycleMatrix(i, :) ...
                                - liftedPathMatrix(incidentVertices(2), :) ...                                
                                - liftedPathMatrix(obj.involutionVertex(incidentVertices(1)), :) ...
                                + liftedPathMatrix(obj.involutionVertex(incidentVertices(2)), :) ...
                                + liftedPathMatrix(incidentVertices(1), :) ;                            
                            cycleMatrix(i, edgesList(j)) = -1;
                            cycleMatrix(i, edgesList(j) + obj.getHalfNumEdges) = 1;
                        end
                        cycleMatrix(i, firstCrossingEdge) = firstCrossingOrientation;
                        cycleMatrix(i, firstCrossingEdge + obj.getHalfNumEdges) = - firstCrossingOrientation;
                    end    
                end    
            end
            prymMatrix = cycleMatrix*transpose(cycleMatrix.*obj.getLength);

        end    

    end
end    