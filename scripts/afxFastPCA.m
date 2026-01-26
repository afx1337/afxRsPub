function [U, explained, mu, V, explainedVoxelPC] = afxFastPCA(y, k, blockSize)
    %[U, explained, mu, V, explainedVoxelPC] = afxFastPCA(y, k, blockSize)
    %
    % Input:
    %   y:                T x v   (Observations x Voxels)
    %   k:                Scalar  Number of components
    %   blockSize:        Scalar  compute covariance/loadings for n voxels at once (default: 1000)
    %
    % Output:
    %   U:                T x k   Scores/Weights
    %   explained         1 x k   Explained variance per component in %
    %   mu:               1 x v   Meain of y
    %   V:                v x k   Spatial Loadings
    %   explainedVoxelPC: v x k   Voxelwise fraction of variance explained by component k
    %
    % Reconstruction:
    %   Yhat = U * V' + mu;
    %
    % comparison to [V,U,~,~,explained,mu] = pca(y,'Economy',true)
    %   - 5-8x faster
    %   - memory efficient
    %   - corr V > 1 - 10e-6
    %   - corr U > 1 - 10e-6
    %   - diff explained < 1.0e-13 (when y is double)
    %   - diff Yhat < 1.0e-08/std (when y is double)

    if nargin < 3 || isempty(blockSize)
        blockSize = 1000;
    end

    [nT, nVox] = size(y);

    % --- 1) Kovarianz im Beobachtungsraum ---
    C = zeros(nT, nT, 'double'); % Kovarianz
    mu = zeros(1,nVox,'like',y); % Mittelwert

    for i = 1:blockSize:nVox
        idx = i:min(i+blockSize-1, nVox);
        mu(idx) = mean(y(:,idx),1,'native');
        Yb = double(y(:,idx) - mu(idx));
        C = C + (Yb * Yb.');
    end

    % --- 2) Eigenzerlegung ---
    [Ud, D] = eigs(C, k, 'largestreal');

    lambda = max(diag(D), eps);
    
    totalVar = trace(C);
    explained = 100 * lambda / totalVar;

    % Scores skalieren (pca()-kompatibel)
    U = cast(Ud .* sqrt(lambda)','like',y);

    % Vorzeichen stabilisieren
    [~,idx] = max(abs(U),[],1);
    signComp = sign(U(sub2ind(size(U),idx,1:k)));
    signComp(signComp==0) = 1;
    U = U .* signComp;

    % --- 3) Spatial Loadings/Varianz berechnen (blockweise) ---
    doV = nargout > 3;
    doExplVox  = nargout > 4;

    if doExplVox
        varVoxel = zeros(nVox,1,'like',y); % Gesamtvarianz je Voxel
        varExpl  = zeros(nVox,k,'like',y); % erklärte Varianz je Komponente
    end
    
    if doV || doExplVox
        invLambda = 1 ./ lambda';
        V = zeros(nVox, k, 'like', y);

        Udbl = double(U);
        for i = 1:blockSize:nVox
            idx = i:min(i+blockSize-1, nVox);
            Yb = double(y(:,idx) - mu(idx));

            % Projection der spatial Loadings
            Vb = (Yb.' * Udbl) .* invLambda;
            V(idx,:) = Vb; % impliziter cast
        
            % voxelweise Varianzmaps
            if doExplVox
                varVoxel(idx) = sum(Yb.^2,1); % Gesamtvarianz je Voxel ohne 1/(N-1), kürzt sich später raus
                varExpl(idx,:) = (Vb.^2) .* lambda'; 
            end
        end
        if doExplVox
            explainedVoxelPC = varExpl ./ max(varVoxel, eps);
        end
    end
end
