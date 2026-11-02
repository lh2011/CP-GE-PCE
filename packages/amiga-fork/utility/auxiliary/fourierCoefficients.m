function [fouierCoefs, dCoefdData] = fourierCoefficients(signal, fftDim)
    L = size(signal, fftDim);
    if ~ismember(fftDim, 1:ndims(signal))
        error("fourierCoefficients: fftDim must in signal size")
    end
    if fftDim ~= 1
        perm = 1:ndims(signal);
        perm(1) = fftDim;
        perm(fftDim) = 1;
        signalReorder = permute(signal, perm);
    else
        signalReorder = signal;
        perm = 1:ndims(signal);
    end
    % check if last values ~= first values
    tol = 1e-6;
    if all(abs(signalReorder(1,:) - signalReorder(end, :)) < tol)
        warning('fourierCoefficients:Periodicity', 'Your first and last values are the same with tolerance %g . For correct Fourier coefficients, do not pass the recurring last value.\n', tol);
    end
    % fft of first dimension always
    mytimer("fft", "start")
    Y = fft(signalReorder, L, 1)/L;
    mytimer("fft", "stop")
    Ncoefs = floor(L/2)+1;
    % truncate signal to coefficients
    fouierCoefs = Y(1:Ncoefs, :);
    % scale to real coefficients
    lastDoubleEntry = ceil(L/2);
    fouierCoefs(2:lastDoubleEntry, :) = 2*fouierCoefs(2:lastDoubleEntry, :);   

    if nargout == 2
        k = (0:Ncoefs-1)';
        m = reshape(0:(L-1), 1, 1, []);
        mytimer("fftder", "start")
        dCoefdData = 1/L.*ones(size(fouierCoefs)).*exp(-1i.*2*pi.*k.*m/L);
        mytimer("fftder", "stop")
        dCoefdData(2:lastDoubleEntry, :, :) = 2*dCoefdData(2:lastDoubleEntry, :, :);
    end

    % create correct sizes again
    sizeY = size(Y);
    fouierCoefs = reshape(fouierCoefs, [size(fouierCoefs,1), sizeY(2:end)]);

    if nargout == 2
        dCoefdData = reshape(dCoefdData, [size(fouierCoefs), L]);
    end

    % permute to original size again
    fouierCoefs = permute(fouierCoefs, perm);

    if nargout == 2
        dCoefdData = permute(dCoefdData, [perm, ndims(signal)+1]); 
    end
end
