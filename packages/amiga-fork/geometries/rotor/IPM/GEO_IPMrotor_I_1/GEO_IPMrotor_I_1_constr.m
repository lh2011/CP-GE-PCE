function c = GEO_IPMrotor_I_1_constr(parameters)

    MAG  = 7e-3;
    MH = 7e-3;
    MW = 19e-3;

    if nargin == 1
        param_names = fieldnames (parameters);
        for iParam  = 1:numel (param_names)
          eval ([param_names{iParam} '= parameters.(param_names{iParam});']);
        end
    end

    c(1) = 3*MW - 2*MAG - 50e-3;
    c(2) = MH + MAG - 15e-3;
    c(3) = MH*MW - 100*1e-6;

    c = c * 100;
end
