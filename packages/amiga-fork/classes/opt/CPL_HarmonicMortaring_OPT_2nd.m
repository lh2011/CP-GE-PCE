classdef CPL_HarmonicMortaring_OPT_2nd < CPL_HarmonicMortaring_OPT & OPT_Element_2nd

    methods
        function Lyy = function_Lyy(obj, t, x, y, z)
            Lyy = sptensor([],[],[obj.NumberDOF, obj.NumberDOF, obj.NumberOptimizationResponses]);
            
            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if ~isempty(indTorque)
                LyyTorque = zeros(obj.NumberDOF, obj.NumberDOF);
                LyyTorque(obj.a2DOFs, obj.couplingDOFs) = obj.Sgn2*obj.MaterialParameter*obj.CouplingInterface2.CouplingMatrixDer *obj.Region2.Length;
                LyyTorque(obj.couplingDOFs, obj.a2DOFs) = obj.Sgn2*obj.MaterialParameter*obj.CouplingInterface2.CouplingMatrixDer'*obj.Region2.Length;
                Lyy(:,:, indTorque) = sptensor(LyyTorque);
            end
        end

        function Lxy = function_Lxy(obj, t, x, y, z)
            Lxy = [];
        end

        function Lxx = function_Lxx(obj, t, x, y, z)
            Lxx = [];
        end
    end
end

