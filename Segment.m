classdef Segment < handle
    % This is a class made to handle subsets of a metric graph which are
    % intervals. All the branching points should be at the boundary of a
    % segment.
    %
    % Concretely a segment is determined by a pair of points that lie on
    % the same edge. The data is stored as a special kind of divisor,
    % but the underlying metric graph is NOT a property of the
    % segment. (This is made to optimize memory usage.)
    %
    % Segments are oriented. The goal of this class is to be used inside
    % specific algorithms requiring to explore a metric graph.
    %
    % See also Divisor, MetricGraph.

    properties (Access = private)
        edgeIndexVector
        distanceVector
        degreeVector
    end
    
    methods (Access = private)
        function addChip(obj, metricGraph, edgeIndex, distanceFromTail, degree)
            if edgeIndex > width(metricGraph.lengths)
                error('The edge index has to be lower than the number of edges of the graph. The given index is %d while the number of edges is %d.',edgeIndex, width(metricGraph.lengths))
            end    
            if edgeIndex < 1
                error('The edge index has to be strictly positive. The given index is %d.', edgeIndex)
            end
            if distanceFromTail < 0
                error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTail)
            end
            if distanceFromTail > metricGraph.lengths(edgeIndex)
                error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTail, metricGraph.lengths(edgeIndex))
            end
            obj.degreeVector = [obj.degreeVector degree];
            obj.distanceVector = [obj.distanceVector [distanceFromTail; metricGraph.lengths(edgeIndex)-distanceFromTail]];
            obj.edgeIndexVector = [obj.edgeIndexVector edgeIndex];
        end

        function reduce(obj, metricGraph)
            if ~isa(metricGraph, 'MetricGraph')
                error('The second input has to be a metric graph, the input class is %s.', class(metricGraph))
            end

       % First step is to take chips that are on a edge to the incident edge of lowest index    
            [obj.edgeIndexVector, obj.distanceVector] = metricGraph.reducePosition(obj.edgeIndexVector, obj.distanceVector);
        % In second, we sort the vector first by edge index then by lengths
        % on the edge
            [~, Sorted] = sort(obj.distanceVector(1, :));
            obj.degreeVector = obj.degreeVector(Sorted);
            obj.distanceVector = obj.distanceVector(:, Sorted);
            obj.edgeIndexVector = obj.edgeIndexVector(Sorted);
            [~, Sorted] = sort(obj.edgeIndexVector);
            obj.degreeVector = obj.degreeVector(Sorted);
            obj.distanceVector = obj.distanceVector(:, Sorted);
            obj.edgeIndexVector = obj.edgeIndexVector(Sorted);
        end
   end

%---------------------------------------------------------------------------------------------    
    
    methods
        % Constructor
        function obj = Segment(metricGraph, edgeIndex, distanceFromTailStart, distanceFromTailEnd)
            % This is the constructor. Segment(metricGraph, edgeIndex,
            % distanceStart, distanceEnd) outputs a segment on the edge
            % edgeIndex of metricGraph which starts at distanceStart and
            % ends at distanceEnd.
            %
            % See also Divisor, MetricGraph.

            if nargin ~=0
                if edgeIndex > width(metricGraph.lengths)
                    error('The edge index has to be lower than the number of edges of the graph. The given index is %d while the number of edges is %d.',edgeIndex, width(metricGraph.lengths))
                end    
                if edgeIndex < 1
                    error('The edge index has to be strictly positive. The given index is %d.', edgeIndex)
                end
                if distanceFromTailStart < 0
                    error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTailStart)
                end
                if distanceFromTailStart > metricGraph.lengths(edgeIndex)
                    error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTailStart, metricGraph.lengths(edgeIndex))
                end
                if distanceFromTailEnd < 0
                    error('The distance from the tail has to be positive. Given distance is %d.', distanceFromTailEnd)
                end
                if distanceFromTailEnd > metricGraph.lengths(edgeIndex)
                    error('The distance from the tail has to be lower than the length of the edge. Given distance is %d while the edge has length %d.', distanceFromTailEnd, metricGraph.lengths(edgeIndex))
                end
            end
            obj.edgeIndexVector = zeros(1,0);
            obj.distanceVector = zeros(2,0);
            obj.degreeVector = zeros(1,0);
            if nargin ~=0
                addChip(obj, metricGraph, edgeIndex, distanceFromTailStart, -1);
                addChip(obj, metricGraph, edgeIndex, distanceFromTailEnd, 1);
            end
        end
                        
        function show(obj)
            % This is the equivalent of the function show of a divisor.
            
            disp([obj.edgeIndexVector; obj.distanceVector; obj.degreeVector])
        end

        function chipVector = getInterior(obj, metricGraph, edgeIndex)
            % mySegment.getInterior(metricGraph, edgeIndex) output the 
            % indices of the points that lies in the interior of the edge 
            % edgeIndex.
            %

            if ~isa(metricGraph, 'MetricGraph')
                error('The input has to be a metric graph, the input class is %s.', class(metricGraph))
            end
            if edgeIndex > width(metricGraph.incidenceMatrix)
                error('You cannot have an edge index larger than the number of edges on the graph. Given index is %d while the number of edges is %d.', edgeIndex, width(metricGraph.incidenceMatrix))
            end
            if edgeIndex < 1 
                error('You cannot have an edge index strictly smaller than 1. Given index is %d.', edgeIndex)
            end
            chipVector = intersect(intersect(find(obj.edgeIndexVector == edgeIndex), find(obj.distanceVector(1, :))), find(obj.distanceVector(2, :)));
        end

        function cloned = clone(obj)
            % mySegment.clone outputs a clone of mySegment.

            cloned = Segment;
            cloned.degreeVector = obj.degreeVector;
            cloned.distanceVector = obj.distanceVector;
            cloned.edgeIndexVector = obj.edgeIndexVector;
        end

        function areSame = areSameUnoriented(obj1, obj2, metricGraph)
            % areSameUnoriented(seg1, seg2, metricGraph) outputs true if
            % seg1 and seg2 have the same ends.

            if ~isa(metricGraph, 'MetricGraph')
                error('The third input has to be a metric graph, the input class is %s.', class(metricGraph))
            end
            Clone1 = clone(obj1);
            Clone1.reduce(metricGraph);
            Clone2 = clone(obj2);
            Clone2.reduce(metricGraph);
            if all(all(Clone1.distanceVector == Clone2.distanceVector)) && all(all(Clone1.edgeIndexVector == Clone2.edgeIndexVector))
                areSame = true;
            else
                areSame = false;
            end    
        end

        function div = turnToDivisor(obj, metricGraph)
            % mySegment.turnToDivisor(metricGraph) outputs a divisor which
            % is -1 at the start of the segment mySegment and 1 at its end.

            if ~isa(metricGraph,'MetricGraph')
                error('The second input has to be a metric graph, the input class is %s.', class(metricGraph))
            end
            div = Divisor(metricGraph);
            if width(obj.edgeIndexVector) ~= 0
                A = find(obj.degreeVector == -1);
                div.addChip(obj.edgeIndexVector(A), obj.distanceVector(1,A), -1)
                A = find(obj.degreeVector == 1);
                div.addChip(obj.edgeIndexVector(A), obj.distanceVector(1,A), 1)
            end    
        end

        function [edge, distance] = getStart(obj)
            % [edge, distance] = getStart(mySegment) will provide the
            % starting point of mySegment.
            %
            % See also getEnd.
  
            startIndex = find(obj.degreeVector == -1);
            edge = obj.edgeIndexVector(startIndex);
            distance = obj.distanceVector(:, startIndex);
        end

        function [edge, distance] = getEnd(obj)
            % [edge, distance] = getEnd(mySegment) will provide the
            % starting point of mySegment.
            %
            % See also getStart.
            endIndex = find(obj.degreeVector == 1);
            edge = obj.edgeIndexVector(endIndex);
            distance = obj.distanceVector(:,endIndex);
        end

        function invert(obj)
            % mySegment.invert will exchange the two extremities of
            % mySegment.

            obj.degreeVector = -obj.degreeVector;
        end

        function towardsHead = isTowardsHead(obj)
            % mySegment.isTowardsHead outputs true if mySegment is oriented
            % in the same direction as the edge it is lying on.

            if obj.edgeIndexVector(1) ~= obj.edgeIndexVector(2)
                error('Segment has to be written on a unique edge')
            end
            if obj.distanceVector(1, find(obj.degreeVector == 1)) > obj.distanceVector(1, find(obj.degreeVector == -1))
                towardsHead = true;
            else
                towardsHead = false;
            end
        end

        function length = lengthToStop(obj, edgeIndex, distance)
            % mySegment.lengthToStop outputs the distance between the
            % starting point of the segment and the boundary of the edge it
            % lies on, where the boundary chosen is on the same side as the
            % end of the segment.
            %
            % mySegment.lengthToStop(edgeIndex, distance) outputs the same
            % data, but the segment is eventually limited by the point
            % described by the pair (edgeIndex, distance).
            %
            % See also changeLengthTo

            if nargin ~= 3
                edgeIndex = 0;
                distance = [0 0];
            end    
            if obj.isTowardsHead
                if obj.edgeIndexVector(1) == edgeIndex && distance(1) >= obj.distanceVector(2, find(obj.degreeVector == -1))
                    length = obj.distanceVector(2, find(obj.degreeVector == -1)) - distance(2);
                else    
                    length = obj.distanceVector(2, find(obj.degreeVector == -1));
                end
            else
                if obj.edgeIndexVector(1) == edgeIndex && distance(1) <= obj.distanceVector(1, find(obj.degreeVector == -1))
                    length = obj.distanceVector(1, find(obj.degreeVector == -1)) - distance(1);
                else    
                    length = obj.distanceVector(1, find(obj.degreeVector == -1));
                end
            end
        end

        function changeLengthTo(obj, length)
            % mySegment.changeLengthTo(length) moves the end point of
            % mySegment sucht that it has a distance length between its two
            % extremities.
            %
            % See also lengthToStop.

            if obj.isTowardsHead
                if obj.distanceVector(2,find(obj.degreeVector == -1)) < length
                    error('Length cannot exceed the distance between the starting point and the end of the edge. Given length is %d while the distance is %d.', length, obj.distanceVector(2,find(obj.degreeVector == -1)))
                end    
                obj.distanceVector(1,find(obj.degreeVector == 1)) = obj.distanceVector(1,find(obj.degreeVector == -1)) + length;
                obj.distanceVector(2,find(obj.degreeVector == 1)) = obj.distanceVector(2,find(obj.degreeVector == -1)) - length;
            else
                if obj.distanceVector(1,find(obj.degreeVector == -1)) < length
                    error('Length cannot exceed the distance between the starting point and the end of the edge. Given length is %d while the distance is %d.', length, obj.distanceVector(1,find(obj.degreeVector == -1)))
                end    
                obj.distanceVector(1,find(obj.degreeVector == 1)) = obj.distanceVector(1,find(obj.degreeVector == -1)) - length;
                obj.distanceVector(2,find(obj.degreeVector == 1)) = obj.distanceVector(2,find(obj.degreeVector == -1)) + length;
            end
        end

        function output = getEdge(obj)
            % mySegment.getEdge outputs the edge the segment lies on.

            if obj.edgeIndexVector(1) == obj.edgeIndexVector(2)
                output = obj.edgeIndexVector(1);
            else
                error('The segment is not represented as on an unique edge, first edge index is %d and second edge index is %d.', obj.edgeIndexVector(1), obj.edgeIndexVector(2))
            end
        end

        function length = getLength(obj)
            % mySegment.getLength outputs the distance between the two
            % extremities of mySegment.
            
            length = abs(obj.distanceVector(1,1) - obj.distanceVector(1,2));
        end    
    end
end