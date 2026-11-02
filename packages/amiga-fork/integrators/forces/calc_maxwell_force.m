function T = calc_maxwell_force(interface, args)
    arguments
        interface {mustBeA(interface, 'IGA_CouplingInterface')};
        args.MaterialParameter = 1/(4*pi*1e-7);
    end

    T = 0;
    a = interface.Region.Solution.Values;
    
    for iBnd = 1:numel(interface.PatchBoundaries)
        patchDOFs = interface.Region.Spaces.gnum{interface.Region.PatchBoundaries(interface.PatchBoundaries(iBnd)).patches};
        A_z = a(patchDOFs,:);
        T = T + args.MaterialParameter * interface.Region.Length * ...
            int_Br_Bt_dGamma(interface.SpacesEval{iBnd}, interface.MeshesEval{iBnd}, A_z);
    end
end