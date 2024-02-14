classdef MetricGraph < handle
    % This is a class made to handle metric graphs in MATLAB.
    %
    % We store our metric graph as a directed graph registered via an 
    % incidence matrix.
    %
    % Lengths are stored as a length vector.
    %
    % The i-th entry of the length vector corresponds to the i-th edge.
    %
    % To input a precise point of the metric graph, the description used is
    % giving first the edge the point lies on and then the distance from
    % the start of the edge to the point. Often the class will prefer
    % receiving a distance vector with first entry being the distance from
    % the start to the point and second entry being the distance from the
    % point to the end of thed edge.

    properties (Access = private)
        numVertices
    end

    properties (Access = protected)
        incidenceMatrix
        lengths
        computedConnected
        isKnownConnected
        computedMultigraph
        isKnownMultigraph
        computedLaplacianCoordinatesX
        computedLaplacianCoordinatesY
        isKnownLaplacianCoordinates
    end

    methods
        % Constructor
        function obj = MetricGraph(numVertices)
            % This is the constructor of MetricGraph.
            % A = MetricGraph(n) creates a metric graph A with n vertices.

            if floor(numVertices) ~= numVertices
                error('The number of vertices has to be an integer. The given number is %d.', numVertices)
            end 
            if numVertices <= 0
                error('The number of vertices has to be at least 1. The given number is %d.',numVertices)
            end
            obj.isKnownConnected = true;
            obj.computedConnected = false;
            obj.isKnownMultigraph = true;
            obj.computedMultigraph = false;
            obj.isKnownLaplacianCoordinates = false;
            obj.numVertices = numVertices;
            obj.incidenceMatrix = zeros(numVertices,0);
            obj.lengths = zeros(1,0);
        end
        
        function addEdge(obj, fromVertex, toVertex, length)
            % myGraph.addEdge(T, H, L) creates an edge to the graph myGraph
            % from the vertex indexed T to the vertex indexed H. 
            % 
            % This edge will have a length L. If no length is precised the
            % edge will have a default length of 1.
            %
            % Edges are ordered by FIFO principle and for practical purposes
            % they are stored as directed edges.
            %
            % See also MetricGraph, removeEdge.
            
            if nargin == 3
                length = 1;
            end    
            if fromVertex > obj.numVertices
                error('The starting vertex number cannot exceed the number of vertices. The given index is %d.',fromVertex)
            end  
            if fromVertex < 1
                error('The starting vertex number cannot be smaller than 1.The given index is %d.',fromVertex)
            end
            if toVertex > obj.numVertices
                error('The ending vertex number cannot exceed the number of vertices. The given index is %d.',toVertex)
            end  
            if toVertex < 1
                error('The ending vertex number cannot be smaller than 1.The given index is %d.',toVertex)
            end
            if fromVertex == toVertex
                error('Self-loops are not allowed.The given indices are %d and %d.',fromVertex, toVertex)
            end
            if length <=0
                error('Length has to be strictly positive. The given length is %d.', length)
            end
            if ~obj.computedConnected
                obj.isKnownConnected = false;
            end    
            if ~obj.computedMultigraph
                obj.isKnownConnected = false;
            end
            obj.isKnownLaplacianCoordinates = false;
            temp = zeros(obj.numVertices,1);
            temp(fromVertex) = -1;
            temp(toVertex) = 1;
            obj.incidenceMatrix = [obj.incidenceMatrix temp];
            obj.lengths = [obj.lengths length];
        end
        
        function length = getLength(obj, numEdge)
            % myGraph.getLength(n) will output the length of the n-th edge
            % of the graph myGraph.
            %
            % myGraph.getLength will output all the lengths.

            if nargin == 1
                length = obj.lengths;
            else
                if numEdge > obj.getNumEdges
                    error('The edge number cannot be larger than the number of edges. The given number is %d.', numEdge)
                end  
                if numEdge < 1
                   error('The edge number cannot be smaller than 1. The given number is %d.', numEdge)
                end 
                length = obj.lengths(numEdge);
            end    
        end

        function incidents = getIncidents(obj, vertex)
            % myGraph.getIncidents(n) will output the indices of the edges
            % that are incident to the vertex numbered n.
            %
            % See also getTailIncidents, getHeadIncidents, getNeighbors.

            if vertex > obj.numVertices
                error('The vertex number cannot exceed the number of vertices. The given index is %d.', vertex)
            end  
            if vertex < 1
                error('The vertex number cannot be smaller than 1. The given index is %d.', vertex)
            end
            incidents = find(obj.incidenceMatrix(vertex, :));
        end

        function incidents = getTailIncidents(obj, vertex)
            % myGraph.getTailIncidents(n) will output the indices of the 
            % edges that start from the vertex number n.
            %
            % See also getIncidents, getHeadIncidents, getNeighbors.

            if vertex > obj.numVertices
                error('The vertex number cannot exceed the number of vertices. The given index is %d.', vertex)
            end  
            if vertex < 1
                error('The vertex number cannot be smaller than 1. The given index is %d.', vertex)
            end
            incidents = find(obj.incidenceMatrix(vertex, :) == -1);
        end

        function incidents = getHeadIncidents(obj, vertex)
            % myGraph.getHeadIncidents(n) will output the indices of the 
            % edges that go to the vertex number n.
            %
            % See also getIncidents, getTailIncidents, getNeighbors.

            if vertex > obj.numVertices
                error('The vertex number cannot exceed the number of vertices. The given index is %d.', vertex)
            end  
            if vertex < 1
                error('The vertex number cannot be smaller than 1. The given index is %d.', vertex)
            end
            incidents = find(obj.incidenceMatrix(vertex, :) == 1);
        end

        function neighbors = getNeighbors(obj, vertex)
            % myGraph.getNeighbors(vertexIndex) will output all the
            % vertices that are adjacent to the vertex vertexIndex. In
            % other words, it will output all the vertices such that there
            % is an edge between these vertices and the vertex vertexIndex.
            %
            % See also getOtherVertex, getIncidents.

            if vertex > obj.numVertices
                error('The vertex number cannot exceed the number of vertices. The given index is %d.', vertex)
            end  
            if vertex < 1
                error('The vertex number cannot be smaller than 1. The given index is %d.', vertex)
            end
            incidentList = getIncidents(obj, vertex);
            tempMatrix = obj.incidenceMatrix(:, incidentList);
            tempMatrix(vertex,:) = 0;
            [neighbors,~] = find(tempMatrix);
        end
       
        function degree = getDegree(obj, vertex)
            % myGraph.getDegree(vertexIndex) outputs the degree of the
            % vertex with index vertexIndex.
            %

            degree = height(obj.getNeighbors(vertex));
        end    


        function numVertices = getNumVertices(obj)
            % myGraph.getNumVertices outputs the number of vertices of the
            % graph myGraph.
            %
            % See also getNumEdges.

            numVertices = obj.numVertices;
        end
      
        function numEdges = getNumEdges(obj)
            % myGraph.getNumEdges outputs the number of edges of the graph 
            % myGraph.
            %
            % See also getNumVertices.
            
            numEdges = width(obj.lengths);
        end

        function tailNumber = getTail(obj, edgeNumber)
            % myGraph.getTail(n) outputs the number of the vertex which is
            % at the start of the n-th edge.
            %
            % See also getHead.

           if edgeNumber > obj.getNumEdges
                error('The edge number cannot be larger than the number of edges. The given number is %d.', edgeNumber)
           end  
           if edgeNumber < 1
                error('The edge number cannot be smaller than 1. The given number is %d.', edgeNumber)
           end    
           tailNumber = find(obj.incidenceMatrix(:,edgeNumber)==-1);
        end
        
        function headNumber = getHead(obj,edgeNumber)
            % myGraph.getHead(n) outputs the number of the vertex which is
            % at the end of the n-th edge.
            %
            % See also getTail.

           if edgeNumber > obj.getNumEdges
                error('The edge number cannot be larger than the number of edges. The given number is %d.', edgeNumber)
           end  
           if edgeNumber < 1
                error('The edge number cannot be smaller than 1. The given number is %d.', edgeNumber)
           end 
           headNumber = find(obj.incidenceMatrix(:,edgeNumber)==1);
        end

        function removeEdge(obj, edgeNumber)
            % myGraph.removeEdge(n) removes the n-th edge of the graph
            % myGraph.
            %
            % See also addEdge, removeAllEdges.

            if edgeNumber > obj.getNumEdges
                error('The edge number cannot be larger than the number of edges. The given number is %d while the number of edges is %d.', edgeNumber, obj.getNumEdges)
            end  
            if edgeNumber < 1
                error('The edge number cannot be smaller than 1. The given number is %d.', edgeNumber)
            end
            if obj.computedConnected
                obj.isKnownConnected = false;
            end
            if obj.computedMultigraph
                obj.isKnownMultigraph = false;
            end
            obj.isKnownLaplacianCoordinates = false;
            obj.lengths(edgeNumber)=[];
            obj.incidenceMatrix(:,edgeNumber)=[];
        end
        
        function removeAllEdges(obj, fromVertex, toVertex)
            % myGraph.removeAllEdges(n,m) removes all the edges that go from
            % the n-th vertex to the m-th vertex.
            %
            % See also addEdge, removeEdge, removeAllEdgesSym.

            if fromVertex > obj.numVertices
                error('The starting vertex number cannot exceed the number of vertices. The given index is %d while the number of vertices is %d.', fromVertex, obj.getNumVertices)
            end  
            if fromVertex < 1
                error('The starting vertex number cannot be smaller than 1.The given index is %d.',fromVertex)
            end
            if toVertex > obj.numVertices
                error('The ending vertex number cannot exceed the number of vertices. The given index is %d while the number of vertices is %d.', toVertex, obj.getNumVertices)
            end  
            if toVertex < 1
                error('The ending vertex number cannot be smaller than 1.The given index is %d.', toVertex)
            end
            if obj.computedConnected
                obj.isKnownConnected = false;
            end
            if obj.computedMultigraph
                obj.isKnownMultigraph = false;
            end
            obj.isKnownLaplacianCoordinates = false;
            obj.removeEdge(find(obj.incidenceMatrix(toVertex,find(obj.incidenceMatrix(fromVertex,:)==-1))==1))
        end

        function removeAllEdgesSym(obj, Vertex1, Vertex2)
            % myGraph.removeAllEdgesSym(n,m) removes all the edges that go
            % between the n-th vertex and the m-th vertex.
            %
            % See also addEdge, removeEdge, removeAllEdges.
            
            if Vertex1 > obj.numVertices
                error('The first vertex number cannot exceed the number of vertices. The given index is %dwhile the number of vertices is %d.', Vertex1, obj.numVertices)
            end  
            if Vertex1 < 1
                error('The first vertex number cannot be smaller than 1.The given index is %d.', Vertex1)
            end
            if Vertex2 > obj.numVertices
                error('The second vertex number cannot exceed the number of vertices. The given index is %d while the number of vertices is %d.', Vertex2, obj.numVertices)
            end  
            if Vertex2 < 1
                error('The second vertex number cannot be smaller than 1.The given index is %d.',fromVertex)
            end
            if obj.computedConnected
                obj.isKnownConnected = false;
            end
            if obj.computedMultigraph
                obj.isKnownMultigraph = false;
            end
            obj.isKnownLaplacianCoordinates = false;
            obj.removeAllEdges(Vertex1, Vertex2)
            obj.removeAllEdges(Vertex2, Vertex1)
        end

        function multigraphic = isMultigraph(obj)
            % myGraph.isMultigraphic returns true if myGraph is a
            % multigraph, i.e. there exists at least one pair of vertices
            % with multiple edges between them.

            if ~obj.isKnownMultigraph
                obj.isKnownMultigraph = true;
                obj.computedMultigraph = false;
                for i = 1:obj.numVertices
                    if  length(obj.getIncidents(i)) ~= length(unique(obj.getNeighbors(i)))
                        obj.computedMultigraph = true;
                    end    
                end
            end
            multigraphic = obj.computedMultigraph;
        end

        function adjMatrix = getAdjacencyMatrixDir(obj)
            % myGraph.getAdjacencyMatrixDir outputs the (asymmetric)
            % adjacency matrix of the metric graph myGraph.
            %
            % See also getAdjacencyMatrix, getLaplacianMatrix, 
            % getIncidenceMatrix.

            adjMatrix = zeros(obj.numVertices);
            for i = 1:obj.getNumEdges
                tail = find(obj.incidenceMatrix(:,i) == -1);
                head = find(obj.incidenceMatrix(:,i) == 1);
                adjMatrix(tail, head) = adjMatrix(tail, head) + 1;
            end
        end
   
        function adjMatrix = getAdjacencyMatrix(obj)
            % myGraph.getAdjacencyMatrixDir outputs the (symmetric)
            % adjacency matrix of the metric graph myGraph.
            %
            % See also getAdjacencyMatrixDir, getLaplacianMatrix,
            % getIncidenceMatrix.
            
            adjMatrix = getAdjacencyMatrixDir(obj);
            adjMatrix = adjMatrix + adjMatrix';
        end
      
        function [edgeIndex, distance] = reducePosition(obj, edgeIndex, distance)
            % [edge, distance] = myGraph.reducePosition(edge, distance)
            % takes as input one point (i.e. an edge and a distance from both
            % ends of the edge, first one being from the tail of the edge
            % and second one from the head of the edge) and outputs the
            % data describing the same point but with an edge index which
            % is minimal.
            %
            % See also vertexToPoint.
            
            if height(distance) ~= 2
                error('The distance vector needs a height of 2. Given vector has height %d.', height(distance))
            end
            if width(edgeIndex) ~= width(distance)
                error('The width of the edge index vector has to be the same as the distance vector. Given width of the first one is %d and the second one is %d.', width(edgeIndex), width(distance))
            end
            if edgeIndex > obj.getNumEdges
                error('The edge index has to be lower than the number of edges. Given number is %d while there are %d edges on the graph.', edgeIndex, obj.getNumEdges)
            end
            if edgeIndex < 1
                error('The edge index has to be strictly positive. Given number is %d.', edgeIndex)
            end
            limitChips = zeros(1,width(distance)) == distance;
            tailSearch = find(limitChips(1,:));
            for i = 1:length(tailSearch)
                pivotVertex = getTail(obj, edgeIndex(tailSearch(i)));
                minIncident = min(getIncidents(obj, pivotVertex));
                if minIncident ~= edgeIndex(tailSearch(i))
                    edgeIndex(tailSearch(i)) = minIncident;
                     if obj.incidenceMatrix(pivotVertex, minIncident) == 1
                         distance(:,tailSearch(i))=[obj.lengths(minIncident);0];
                     else
                         distance(:,tailSearch(i))=[0;obj.lengths(minIncident)];
                     end    
                end    
             end

             headSearch = find(limitChips(2,:));
             for i = 1:length(headSearch)
                pivotVertex = getHead(obj, edgeIndex(headSearch(i)));
                minIncident = min(getIncidents(obj, pivotVertex));
                if minIncident ~= edgeIndex(headSearch(i))
                   edgeIndex(headSearch(i)) = minIncident;
                   if obj.incidenceMatrix(pivotVertex, minIncident) == 1
                       distance(:,headSearch(i))=[obj.lengths(minIncident);0];
                   else
                       distance(:,headSearch(i))=[0;obj.lengths(minIncident)];
                   end    
                end    
             end
        end

       function cloned = clone(obj)
           % myGraph.clone outputs a copy of the graph myGraph, which will
           % have exactly the same properties.

           if ~isa(obj,'MetricGraph')
                error('The input has to be a metric graph, the input class is %s.', class(obj))
            end    
            cloned = MetricGraph(obj.numVertices);
            cloned.incidenceMatrix = obj.incidenceMatrix;
            cloned.lengths = obj.lengths;
            cloned.isKnownConnected = obj.isKnownConnected;
            cloned.isKnownMultigraph = obj.isKnownMultigraph;
            cloned.computedConnected = obj.computedConnected;
            cloned.computedMultigraph = obj.computedMultigraph;
            cloned.isKnownLaplacianCoordinates = obj.isKnownLaplacianCoordinates;
            cloned.computedLaplacianCoordinatesX = obj.computedLaplacianCoordinatesX;
            cloned.computedLaplacianCoordinatesY = obj.computedLaplacianCoordinatesY;
       end

       function plot(obj, X, Y)
           % myGraph.plot(X, Y) creates a visual plot of myGraph where the
           % i-th vertex is placed at coordinates (X(i), Y(i)).
           %
           % myGraph.plot executes myGraph.plot(X, Y) with X and Y
           % being the laplacian coordinates.
           %
           % See also tikzPlot, getLaplacianCoordinates.
        
            if nargin == 1
                [X,Y] = obj.getLaplacianCoordinates;
            end
            if obj.numVertices ~= width(X)
                error('The width of the vector X has to be the number of vertices of the graph. The graph has %d vertices and the vector has width %d.',obj.numVertices, width(X))
            end
            if obj.numVertices ~= width(Y)
                error('The width of the vector Y has to be the number of vertices of the graph. The graph has %d vertices and the vector has width %d.',obj.numVertices, width(Y))
            end
            % Create an empty plot
            clf
            % figure;
            hold on;
            % Plot edges
            for i = 1:obj.getNumEdges
                fromVertex = find(obj.incidenceMatrix(:,i) == -1);
                toVertex = find(obj.incidenceMatrix(:,i) == 1);
                plot([X(fromVertex), X(toVertex)], [Y(fromVertex), Y(toVertex)], 'k', 'LineWidth', 2, 'Marker','.')
            end
            % Plot each vertex as a node
            for vertex = 1:obj.numVertices
                plot(X(vertex), Y(vertex), '.', 'MarkerSize', 20, 'Color','black');
            end             
            % Set axis properties
            axis([min(X)-1 max(X)+1 min(Y)-1 max(Y)+1]);
            axis equal;
            axis off;
            hold off;
       end

       function tikzPlot(obj, X, Y, needLastTikzLine)
           % myGraph.tikzPlot(X, Y, needLastTikzLine) creates the Tikz code 
           % that produces the visual plot of the graph myGraph where the
           % i-th vertex is located at coordinates (X(i), Y(i)).
           %
           % myGraph.tikzPlot executes myGraph.tikzPlot(X, Y) with X and Y
           % being the laplacian coordinates.
           % 
           % One can ommit the input needLastTikzLine. If this input is set 
           % to false, the code doesn't produce the \end{tikzpicture}.
           %
           % The Tikz code is put in the file named "output.txt", in your
           % current folder. If the file doesn't exist then it is created.
           %
           % See also plot, getLaplacianCoordinates.
           
            if nargin == 1
                [X,Y] = obj.getLaplacianCoordinates;
            end
            if obj.numVertices ~= width(X)
                error('The width of the vector X has to be the number of vertices of the graph. The graph has %d vertices and the vector has width %d.',obj.numVertices, width(X))
            end
            if obj.numVertices ~= width(Y)
                error('The width of the vector Y has to be the number of vertices of the graph. The graph has %d vertices and the vector has width %d.',obj.numVertices, width(Y))
            end    
            if nargin <= 3
                needLastTikzLine = true;
            end
            if ~isa(needLastTikzLine, 'logical')
                error('The fourth input has either to be empty or logical. The input class is %s.', class(needLastTikzLine))
            end    
            fileID = fopen("output.txt",'a+');
            fprintf(fileID, '%%--------------------------\n');
            fprintf(fileID, '\\begin{tikzpicture}\n');
            % Plot edges
            for i = 1:obj.getNumEdges
                fromVertex = find(obj.incidenceMatrix(:,i) == -1);
                toVertex = find(obj.incidenceMatrix(:,i) == 1);
                fprintf(fileID, '    \\draw[black] (%d, %d) -- (%d, %d); \n', X(fromVertex),Y(fromVertex), X(toVertex),Y(toVertex));
            end
             % Plot each vertex as a node
            for vertex = 1:obj.numVertices
                fprintf(fileID, '    \\filldraw[black, thick] (%d, %d) circle (2pt); \n', X(vertex),Y(vertex));
            end
            if needLastTikzLine
                fprintf(fileID, '\\end{tikzpicture}\n');
            end    
            fclose(fileID);
        end

        function otherVertex = getOtherVertex(obj, edgeIndex, vertex)
            % myGraph.getOtherVertex(edge, vertexIndex) outputs the index
            % of the vertex which is incident to the given edge but which
            % is not the vertex numbered vertexIndex. If the given vertex
            % is not incident, both incident vertices will be outputed.
            %
            % See also getNeighbors, getIncidentVertices.
            
            if edgeIndex > obj.getNumEdges
                error('The edge index has to be lower than the number of edges. Given number is %d while there are %d edges on the graph.', edgeIndex, obj.getNumEdges)
            end
            if edgeIndex < 1
                error('The edge index has to be strictly positive. Given number is %d.', edgeIndex)
            end
            if vertex > obj.numVertices
                error('The vertex number cannot exceed the number of vertices. The given index is %d.', vertex)
            end  
            if vertex < 1
                error('The vertex number cannot be smaller than 1. The given index is %d.', vertex)
            end
            [otherVertex, ~] = find(obj.incidenceMatrix(:, edgeIndex));
            otherVertex(otherVertex == vertex) = [];
        end

        function otherVertex = getIncidentVertices(obj, edgeIndex)
            % myGraph.getIncidentVertices(edgeIndex) outputs the index
            % of the vertex which are incident to the given edge.
            %
            % See also getNeighbors, getOtherVertex.
            
            if edgeIndex > obj.getNumEdges
                error('The edge index has to be lower than the number of edges. Given number is %d while there are %d edges on the graph.', edgeIndex, obj.getNumEdges)
            end
            if edgeIndex < 1
                error('The edge index has to be strictly positive. Given number is %d.', edgeIndex)
            end            
            [otherVertex, ~] = find(obj.incidenceMatrix(:, edgeIndex));
        end

        function [edge, distance] = vertexToPoint(obj, vertex)
            % [edge, distance] = myGraph.vertexToPoint(n) outputs a pair
            % edge-distance describing the n-th vertex of the graph
            % myGraph.
            

            edge = obj.getIncidents(vertex);
            edge = edge(1);
            if obj.incidenceMatrix(vertex,edge) == -1
                distance = [0; obj.lengths(edge)];
            else
                distance = [obj.lengths(edge); 0];
            end            
        end

        function connected = isConnected(obj)
            % myGraph.isConnected outputs true if the metric graph myGraph 
            % is a (path) connected topological space.
            %

            if ~obj.isKnownConnected                         
                obj.isKnownConnected = true;
                toExplore = obj.getIncidents(1);
                toExploreOrigins = ones(1,width(toExplore));
                explored = toExplore;
                listOfReachedVertices = zeros(obj.getNumVertices,1);
                listOfReachedVertices(1) = true;
                while width(toExplore) ~= 0
                    currentVertex = obj.getOtherVertex(toExplore(1), toExploreOrigins(1));
                    listOfReachedVertices(currentVertex) = 1;
                    listOfAdjacents = obj.getIncidents(currentVertex);
                    listOfAdjacents = setdiff(listOfAdjacents, intersect(listOfAdjacents, explored));
                    toExplore = [toExplore listOfAdjacents]; %#ok
                    explored = [explored listOfAdjacents]; %#ok
                    toExploreOrigins = [toExploreOrigins currentVertex*ones(1,width(listOfAdjacents))]; %#ok                             
                    toExplore(1) = [];
                    toExploreOrigins(1) = [];
                    if all(listOfReachedVertices)
                        break
                    end    
                end
                obj.computedConnected = all(listOfReachedVertices);
            end
            connected = obj.computedConnected;
        end

        function changeOrientation(obj, edgeIndex)
            % myGraph.changeOrientation(edgeIndex) changes the orientation
            % of the edge numbered edgeIndex.
            %

            if edgeIndex > obj.getNumEdges
                error('The edge number cannot be larger than the number of edges. The given number is %d.', edgeNumber)
            end  
            if edgeIndex < 1
                error('The edge number cannot be smaller than 1. The given number is %d.', edgeNumber)
            end
            obj.incidenceMatrix(:,edgeIndex) = -obj.incidenceMatrix(:,edgeIndex);
        end

        function matrix = getLaplacianMatrix(obj)
            % myGraph.getLaplacianMatrix outputs the laplacian matrix of
            % the graph myGraph.
            %
            % See also getAdjacencyMatrix, getAdjacencyMatrixDir,
            % getLaplacianCoordinates, getIncidenceMatrix.

            matrix = -obj.getAdjacencyMatrix;
            for i = 1:obj.getNumVertices
                matrix(i,i) = obj.getDegree(i);
            end    
        end

        function matrix = getIncidenceMatrix(obj, row, col)
            % myGraph.getIncidenceMatrix outputs the incidence matrix of
            % the graph myGraph. Column number n contains -1 at the index
            % of the vertex where the n-th edge starts and 1 at the index
            % of the vertex where the n-th edge ends.
            %
            % myGraph.getIncidenceMatrix(A, B) the entry (A, B) of the
            % incidence matrix.
            %
            % See also getAdjacencyMatrix, getAdjacencyMatrixDir,
            % getLaplacianMatrix.

            if nargin == 3
                matrix = obj.incidenceMatrix(row, col);
            else
                matrix = obj.incidenceMatrix;
            end    
        end

        function [X, Y] = getLaplacianCoordinates(obj, scaleFactor)
            % [X, Y] = myGraph.getLaplacianCoordinates(scaleFactor) outputs
            % the second and third eigenvectors of the laplacian matrix of
            % myGraph. These vectors are normalized such that the largest
            % number appearing has absolute value of scaleFactor.
            %
            % If scaleFactor is not precised, the default value used is 3.
            %
            % See also getLaplacianMatrix, plot, tikzPlot.
            
            if nargin < 2
                scaleFactor = 3;
            end
            if scaleFactor <= 0
                error('Scalefactor has to be strictly positive, given input is %d.', scaleFactor)
            end
            if ~obj.isKnownLaplacianCoordinates
                obj.isKnownLaplacianCoordinates = true;
                switch obj.getNumVertices
                    case 1
                        obj.computedLaplacianCoordinatesX = 0;
                        obj.computedLaplacianCoordinatesY = 0;
                    case 2
                        obj.computedLaplacianCoordinatesY = [0 0];
                        obj.computedLaplacianCoordinatesX = [-1 1];
                    otherwise
                        [eigenVectors, eigenValues] = eig(-obj.getLaplacianMatrix);
                        [~, order] = sort(diag(eigenValues), 'descend');
                        obj.computedLaplacianCoordinatesX = transpose(eigenVectors(:, order(2)));
                        obj.computedLaplacianCoordinatesY = transpose(eigenVectors(:, order(3)));
                end              
            end
            scale = max(max(abs(obj.computedLaplacianCoordinatesX)), max(abs(obj.computedLaplacianCoordinatesY)));
            if scale ~= 0
                X = scaleFactor*obj.computedLaplacianCoordinatesX/scale;
                Y = scaleFactor*obj.computedLaplacianCoordinatesY/scale;
            else
                X = obj.computedLaplacianCoordinatesX;
                Y = obj.computedLaplacianCoordinatesY;
            end
        end

        function canDiv = canonicalDivisor(obj)
            % myGraph.canonicalDivisor outputs the canonical divisor of
            % myGraph.
            
            if ~obj.isConnected
                error('The input metric graph has to be connected.')
            end
            canDiv = Divisor(obj);
            for i = 1:obj.getNumVertices
                [edge, distance] = obj.vertexToPoint(i);
                canDiv.addChip(edge, distance, obj.getDegree(i) - 2)
            end    
        end

        function [cycleMatrix, supportingTree, pathMatrix] = getBasisOfCycles(obj)
            % [cycleMatrix, supportingTree, pathMatrix] =
            % myGraph.getBasisOfCycles outputs three informations. 
            % 
            % The cycleMatrix will contain a basis of the indenpendent cycles
            % of myGraph. Every row will correspond to a cycle, every
            % column to an edge and a +1 tells that the edge is present in
            % the cycle in its orientation and a -1 that the edge is
            % present but with its opposite orientation.
            %
            % The supportingTree is a row vector containing 1 at the
            % indices of the supporting tree used to find the basis of
            % cycles.
            %
            % The pathMatrix is a matrix such that the n-th row contains the
            % path inside the supporting tree from the vertex 1 to the
            % n-th vertex.
            

            if ~obj.isConnected
                error('The input metric graph has to be connected.')
            end

            pointer = 1;
            pathMatrix = zeros(obj.getNumVertices, obj.getNumEdges);
            supportingTree = zeros(1, obj.getNumEdges);
            cycleMatrix = zeros(obj.getNumEdges - obj.getNumVertices + 1, obj.getNumEdges);
            toExplore = obj.getIncidents(1);
            toExploreOrigins = ones(1,width(toExplore));
            explored = toExplore;
            while width(toExplore) ~= 0
                currentVertex = obj.getOtherVertex(toExplore(1), toExploreOrigins(1));
                listOfAdjacents = obj.getIncidents(currentVertex);
                listOfAdjacents = setdiff(listOfAdjacents, intersect(listOfAdjacents, explored));
                toExplore = [toExplore listOfAdjacents]; %#ok
                explored = [explored listOfAdjacents]; %#ok
                toExploreOrigins = [toExploreOrigins currentVertex*ones(1,width(listOfAdjacents))]; %#ok               
                % The next if will add the information about the vertex to
                % our matrices containing the paths.
                if any(pathMatrix(currentVertex,:))
                    loopExplored = pathMatrix(toExploreOrigins(1),:) - pathMatrix(currentVertex,:);
                    if obj.getTail(toExplore(1)) == toExploreOrigins(1)
                        loopExplored(toExplore(1)) = 1;
                    else
                        loopExplored(toExplore(1)) = -1;
                    end
                    cycleMatrix(pointer,:) = loopExplored;
                    pointer = pointer + 1;
                else
                    supportingTree(toExplore(1)) = 1;
                    pathMatrix(currentVertex,:) = pathMatrix(toExploreOrigins(1),:);
                    if obj.getTail(toExplore(1)) == toExploreOrigins(1)
                        pathMatrix(currentVertex, toExplore(1)) = 1;
                    else
                        pathMatrix(currentVertex, toExplore(1)) = -1;
                    end    
                end
                toExplore(1) = [];
                toExploreOrigins(1) = [];
            end
        end

        function genus = getGenus(obj)
            % myGraph.getGenus outputs the genus of myGraph. In other
            % words, it outputs the number of edges minus the number of
            % vertices + 1.
            %
            % The graph needs to be connected.

            if ~obj.isConnected
                error("The input graph is not connected.")
            end
            genus = obj.getNumEdges - obj.getNumVertices + 1;
        end

        function jacobianMatrix = getJacobianVariety(obj)
            % myGraph.getJacobianVariety outputs the Jacobian variety of
            % the graph. It gives the matrix of the bilinear form
            % corresponding to the homology of myGraph.
            %
            % See also getBasisOfCycles.

            [cycleMatrix, ~] = obj.getBasisOfCycles;
            jacobianMatrix = cycleMatrix*transpose(cycleMatrix.*obj.getLength);
        end    

    end
end