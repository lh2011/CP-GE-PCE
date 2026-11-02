classdef EXC_CentrifugalForce < EXC
    properties
        
    end

    methods
        function obj = EXC_CentrifugalForce(material, regionNumbers, omega)
            obj@EXC(material, regionNumbers);
            m_r_omega2 = @(x,y) material.getRho()*omega^2 .* cat(1, reshape (cos(atan2(y,x)).* sqrt(x.^2 + y.^2), [1, size(x)]), reshape (sin(atan2(y,x)).* sqrt(x.^2 + y.^2), [1, size(x)]));
            obj.setValue(m_r_omega2);
        end
    end
end