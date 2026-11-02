

% integration of losses according to the Steinmetz formula
% P = Pe + Ph = int{ ke*(f*B)^2 + kh*f^alpha*B^beta }


% deltaT, time length of signal, can be vectorized
function varargout = int_steinmetz(spu, msh, deltaT, Az, args)
    arguments
        spu
        msh
        deltaT (:,1) double
        Az 
        % https://ieeexplore.ieee.org/document/9944887
        args.ke = 6.36e-5;      % eddy current loss factor [W/kg]
        args.kh = 7.55e-3;      % hysteresis loss factor [W/kg]
        args.alpha = 1.3;       % 
        args.beta = 1.8;        %
        args.density = 7650;    % [kg/m^3]
        args.length = 1;        % axial length [m]
    end

    % Ndof x Nel x nT
    Az_el = zeros([size(spu.connectivity), size(Az, 2)]);
    for iel = 1:msh.nel
        Az_el(:, iel, :) = Az(spu.connectivity(:, iel), :);
    end

    gradu = permute(spu.shape_function_gradients, [1,2,3,5,4]); % [dim, nquad, u, t, nel]
    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);
    gradNi_ai = sum(gradu.*permute(Az_el, [4,5,1,3,2]), 3);
    
    B_mag = (sum(gradNi_ai.^2, 1)).^0.5;
    Bcoeff = fourierCoefficients(B_mag, 4);
    % fouerier coefficients
    Bhat = abs(Bcoeff);
    % ground frequencies
    fs = 1./reshape(deltaT, [], 1);
    % fourier frecquencies
    fourierfs = reshape(0:(size(Bhat, 4)-1), 1, 1, 1, []);
    f = fs.*fourierfs;
    
    Pe = args.density*args.length * sum(jacdet.*args.ke.*(f.*Bhat).^2, [2,3,4,5]);
    Ph = args.density*args.length * sum(jacdet.*args.kh.*f.^args.alpha.*Bhat.^args.beta, [2,3,4,5]);

    if nargout == 1
        varargout{1} = Pe + Ph;
    elseif nargout == 2
        varargout{1} = Pe;
        varargout{2} = Ph;
    else
        error("Wrong number of output arguments");
    end

    % %% plot for validation
    % fMax = size(Bcoeff, 4);
    % T = deltaT(end);
    % Bcoeff = Bcoeff(:,:,:,1:fMax,:);
    % t_test = reshape(linspace(0, T, 200), 1, 1, 1, 1, 1, []);
    % freqs = reshape(0:size(Bcoeff, 4)-1, 1, 1, 1, [])/T;
    % 
    % Brec = sum(real(Bcoeff.*exp(1i.*2*pi*freqs.*t_test)), 4);
    % 
    % figure(123)
    % clf
    % subplot(2,1,1);
    % q = randi(size(Bcoeff, 2));
    % e = randi(size(Bcoeff, 5));
    % time = linspace(0, T, size(B_mag, 4)+1);
    % time = time(1:end-1);
    % plot(time, squeeze(B_mag(1,q,1,:,e)))
    % xlabel("time");
    % ylabel("B amplitude");
    % title(sprintf("B of element %i at quadrature point %i", e, q))
    % hold on
    % plot(squeeze(t_test), squeeze(Brec(1,q,1,1,e,:)))
    % subplot(2,1,2);
    % stem(0:size(Bhat, 4)-1,squeeze(Bhat(:, q, :, :, e)));
    % xlabel("fourier coefficient")
    % ylabel("B")

end
