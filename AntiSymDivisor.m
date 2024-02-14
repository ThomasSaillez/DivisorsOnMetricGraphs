classdef AntiSymDivisor < handle & Divisor
    methods
        function obj = AntiSymDivisor(freeDoubleCover)
            % This is the constructor of the class AntiSymDivisor.
            % D = Divisor(A) creates a divisor defined on the free double 
            % covering A.
            %
            % See also Divisor, FreeDoubleCovering.

            if ~isa(freeDoubleCover,'FreeDoubleCovering')
                error('The input has to be a metric graph, the input class is %s.', class(freeDoubleCover))
            end
            if ~freeDoubleCover.isConnected
                error('The input free double covering has to be connected and it is not.')
            end    
            obj@Divisor(freeDoubleCover);
            obj.edgeIndexVector = zeros(1,0);
            obj.distanceVector = zeros(2,0);
            obj.degreeVector = zeros(1,0);
        end

        function addChip(obj, edgeIndex, distanceFromTail, degree)
            % myDivisor.addChip(edgeIndex, distanceFromTail, degree) will
            % increase the divisor myDivisor by degree at the position
            % (edgeIndex, distanceFromTail) and decrease its value by degree 
            % at the point which is the involution of the previous point.
            % Then it will try to simplify the storage of the divisor.
            %
            % See also removeChip.

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
            [edgeIndex, distanceFromTail] = obj.metricGraph.involution(edgeIndex, distanceFromTail);
            addChipPrivate(obj, edgeIndex, distanceFromTail, -degree);
            reduceVectors(obj);
        end
        
        function removeChip(obj, chipIndex)
            % myDivisor.removeChip(chipIndex) removes the chip of index
            % chipIndex of myDivisor and its antisymmetric image.
            %
            % See also addChip and show.

            if chipIndex > width(obj.degreeVector)
                error('You cannot have a chip index higher than number of chips. Given index is %d and the number of chips is %d.', chipIndex, width(obj.degreeVector))
            end 
            if chipIndex < 1
                error('You cannot have a chip index lower than 1. Given index is %d.', chipIndex)
            end  
            obj.addChip(obj.edgeIndexVector(chipIndex), obj.distanceVector(1, chipIndex), -obj.degreeVector(chipIndex))
        end

        function cloned = clone(obj)
            % myDivisor.clone outputs a copy of the antisymmetric divisor
            % myDivisor.

            if ~isa(obj,'AntiSymDivisor')
                error('The input has to be a AntiSymDivisor, the input class is %s.', class(obj))
            end
            cloned = AntiSymDivisor(obj.metricGraph);
            cloned.degreeVector = obj.degreeVector;
            cloned.distanceVector = obj.distanceVector;
            cloned.edgeIndexVector = obj.edgeIndexVector;
        end

        function divisor = forgetAntiSymmetry(obj)
            % myDivisor.forgetAntiSymmetry outputs myDivisor has a 
            % Divisor.

            divisor = Divisor(obj.metricGraph);
            divisor.degreeVector = obj.degreeVector;
            divisor.distanceVector = obj.distanceVector;
            divisor.edgeIndexVector = obj.edgeIndexVector;
        end

        function degree = antiSymmetricDegree(obj)
            % myDivisor.antiSymmetricDegree outputs the total sum 
            % of the positive degrees of the individual chips
            % of the divisor.
            %
            % See also degreeOfPoint.

            degree = (obj.degreeVector > 0) * transpose (obj.degreeVector);
        end
    end
end    