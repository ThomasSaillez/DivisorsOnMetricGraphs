classdef DoubleCovering < MetricGraph & handle
    % This is a class made to handle double coverings of metric graphs.
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
    properties
        isCrossing
    end
    methods
        %Constructor
        function obj = DoubleCovering(metricGraph)
            % This is the constructor of DoubleCovering.
            % myCover = DoubleCovering(myGraph) makes a double covering out
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
            obj.isCrossing = false(metricGraph.getNumVertices, 1);
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

            metricGraph = MetricGraph(obj.numVertices);
            metricGraph.incidenceMatrix = obj.incidenceMatrix;
            metricGraph.lengths = obj.lengths;
        end

        function [edgeIndex, distance] = involution(obj, edgeIndex, distance)
            % myCover.involution(edgeIndex, distance) computes the image of
            % the point described by the pair (edgeIndex, distance) under 
            % the special involution of myCover.
            
            if edgeIndex > obj.getNumEdges
                error('The edge index has to be lower than the number of edges of the graph. The given index is %d while the number of edges is %d.', edgeIndex, obj.getNumEdges)
            end    
            if edgeIndex < 1
                error('The edge index has to be strictly positive. The given index is %d.', edgeIndex)
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

        function addEdge(obj) %#ok
            % One may not add an edge to a double covering.
            %
            % See also forgetCovering.
            error('One may not add an edge to a double covering.')
        end

        function removeEdge(obj) %#ok
            % One may not remove an edge from a double covering.
            %
            % See also forgetCovering.
            error('One may not remove an edge from a double covering.')
        end

        function removeAllEdges(obj) %#ok
            % One may not remove an edge from a double covering.
            %
            % See also forgetCovering.
            error('One may not remove an edge from a double covering.')
        end

        function removeAllEdgesSym(obj) %#ok
            % One may not remove an edge from a double covering.
            %
            % See also forgetCovering.
            error('One may not remove an edge from a double covering.')
        end
        
        function changeOrientation(obj) %#ok
            % One may not change an edge from a double covering.
            %
            % See also forgetCovering.
            error('One may not change an edge from a double covering.')
        end

        function cloned = clone(obj)
           % myCover.clone outputs a copy of the graph myGraph, which will
           % have exactly the same properties.

           if ~isa(obj,'DoubleCovering')
                error('The input has to be a double covering, the input class is %s.', class(obj))
           end    
           cloned = DoubleCovering(obj.baseGraph);
           cloned.incidenceMatrix = obj.incidenceMatrix;
           cloned.lengths = obj.lengths;
           cloned.isCrossing = obj.isCrossing;
       end
    end
end    