classdef PMSM_MotorClass_OPT_2nd < PMSM_MotorClass_OPT & OPT_System_2nd

    methods
        function obj = PMSM_MotorClass_OPT_2nd(rotor, stator, coupling)
            obj@PMSM_MotorClass_OPT(rotor, stator, coupling);
        end

        function H = hessian_responses(obj, x, t)

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            assert(~isempty(indTorque), "Torque needs to be a response in order to calculate mean and std!");
            indTorqueMean = findStructName(obj.OptimizationResponses, "TorqueMEAN");
            indTorqueStd = findStructName(obj.OptimizationResponses, "TorqueSTD");
            
            Torques = zeros(numel(obj.OptimizationAngles), 1);
            dTorquesdx = zeros(numel(obj.OptimizationAngles), obj.NumberOptimizationParameters);
            dTorquesdxdx = zeros(numel(obj.OptimizationAngles), obj.NumberOptimizationParameters, obj.NumberOptimizationParameters);

            for iAngle = 1:numel(obj.OptimizationAngles)
                theta = obj.OptimizationAngles(iAngle);
                obj.Coupling.setRotationAngle(theta);
                
                [Hi, gradri, ri] = hessian_responses@OPT_System_2nd(obj, x, theta);
                
                Torques(iAngle) = ri(indTorque);
                dTorquesdx(iAngle, :) = gradri(indTorque, :);
                dTorquesdxdx(iAngle, :, :) = Hi(indTorque, :, :);
                % Use values of first rotation angle for all other
                % responses which are not mean or std torque
                if iAngle == 1
                    H = Hi;
                end
            end

            obj.LastTorques = Torques;
            obj.LastTorquesx = x;

            Tmean = mean(Torques);
            Tstd = std(Torques, 1, 1);
            
            N = numel(obj.OptimizationAngles);
            dTdPl = dTorquesdx;
            dTdPm =  permute(dTdPl, [1,3,2]);
            dTmeandPl = mean(dTdPl, 1);
            dTmeandPm = permute(dTmeandPl, [1,3,2]);

            dTstddPl = 1/Tstd*(1/N*sum(Torques.*dTdPl, 1) - Tmean*dTmeandPl);
            dTstddPm = permute(dTstddPl, [1,3,2]);

            dTmeandPldPm = mean(dTorquesdxdx, 1);

            dTstddPldPm = 1/Tstd*(-dTstddPl.*dTstddPm + 1/N*(sum(dTdPm.*dTdPl, 1) ...
                + sum(Torques.*dTorquesdxdx, 1)) - dTmeandPl.*dTmeandPm - Tmean*dTmeandPldPm);

            if indTorqueMean > 0
                H(indTorqueMean, :,:) = squeeze(dTmeandPldPm); 
            end

            if indTorqueStd > 0
                H(indTorqueStd, :,:) = squeeze(dTstddPldPm); 
            end
        end    
    end
end
