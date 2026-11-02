% return list of patches and parametric points for physical point

function [xParam, patchNumbers] = physical2param(nurbs, xPhysical)

    xParam = zeros(size(xPhysical));
    patchNumbers = zeros(1, size(xPhysical, 2));

    % reference map for good initial condition / patch guess
    xrefParam = [0.1,0.5,0.9];
    yrefParam = [0.1,0.5,0.9];
    nPatches = numel(nurbs);
    nRef = numel(xrefParam)*numel(yrefParam);
    xRef = zeros(2, numel(nurbs), nRef);
    for nrbNumber = 1:nPatches  
        xrefPhysical = nrbeval(nurbs(nrbNumber), {xrefParam,yrefParam});
        xRef(:, nrbNumber, :) = reshape(xrefPhysical(1:2, :, :), 2, [], 1);
    end
    progressbar("Inverse Mapping for Points:")
    for iPoint = 1:size(xPhysical, 2)
        progressbar( iPoint/size(xPhysical, 2));
        xPhy = xPhysical(:, iPoint);
        xDiff = reshape(vecnorm(xPhy - xRef, 2, 1), nPatches, nRef);
        % minimum difference to forward values for all patches
        [minNurbs, minIndNurbs] = min(xDiff, [], 2);
        [~, testOrder] = sort(minNurbs);

        % ordering of all nurbs
        for nrbNumber = reshape(testOrder, 1, [])
            [indxp, indyp] = ind2sub([numel(xrefParam), numel(yrefParam)], minIndNurbs(nrbNumber));
            xStartParam = [xrefParam(indxp); yrefParam(indyp)];
            [xPar, converge] = nrbinverse(nurbs(nrbNumber), [xPhy; 0], 'u0', xStartParam, 'Display', false);
            if (converge)
                patchNumbers(iPoint) = nrbNumber;
                xParam(:,iPoint) = xPar;
                break;
            end
            if nrbNumber == testOrder(end)
                warning("physical2param:PointNotFound", "Point Number %i not found", iPoint);
            end
        end
    end
    progressbar("done");
    not_found = find(patchNumbers==0);
    if (not_found)
        xParam = xParam(:, patchNumbers ~= 0);
        patchNumbers = patchNumbers(patchNumbers ~= 0);
        warning('physical2param:PointNotFound', 'Your physical points: %s were not found in the computation domain!', join(string(not_found), " "));
    end
end
