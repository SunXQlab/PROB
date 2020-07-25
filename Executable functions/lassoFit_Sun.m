% ===================================================
%                 lassoFit_Sun() 
% ===================================================
function [B,Intercept,lambda,mspe] = ...
    lassoFit_Sun(X,Y,weights,lambda,alpha,Intercept_Flag,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE,maxIter)
%
% ------------------------------------------------------
% Perform model fit for each lambda and the given alpha
% ------------------------------------------------------

[N,P] = size(X);
nLambda = length(lambda);

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ever_active & ~constantPredictors;

% === standardization and weights ===
%
observationWeights = ~isempty(weights);
if ~isempty(weights)
    observationWeights = true;
    weights = weights(:)';
    % Normalize weights up front.
    weights = weights / sum(weights);
end

if ~observationWeights
    muY = mean(Y);
else
    muY = weights*Y;
end
Y0 = bsxfun(@minus,Y,muY);

if standardize
    if ~observationWeights
        % Center and scale
        [X0,muX,sigmaX] = zscore(X,1);
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
    else
        % Weighted center and scale
        muX = weights*X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = sqrt( weights*(X0.^2) );
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
        X0 = bsxfun(@rdivide, X0, sigmaX);
    end
else
    if ~observationWeights
        % Center
        muX = mean(X,1);
        X0 = bsxfun(@minus,X,muX);
        sigmaX = 1;
    else
        % Weighted center
        muX = weights*X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = 1;
    end
end

% If using observation weights, make a weighted copy of the 
% predictor matrix, to save time in the weighted partial regressions.
if observationWeights
    wX0 = bsxfun(@times, X0, weights');
    totalweight = 1;
else
    wX0 = X0;
    totalweight = N;
end

% b will be the current coefficient estimate, iteratively updated.
% Because we retain b from one value of lambda to the next,
% we get a de facto warm start.
b = zeros(P,1);

% Preallocate the returned matrix of coefficients, B, and the intercepts.
B = zeros(P,nLambda);

active = false(1,P);

for i = 1:nLambda
    
    lam = lambda(i);
    if lam >= lambdaMax
        continue;
    end
    threshold = lam * alpha;
    
    % Denominator in coordinate descent update
    if standardize
        if observationWeights
            shrinkFactor = weights*(X0.^2) + lam*(1 - alpha);
        else
            shrinkFactor = ones(1,P) * (1 + lam*(1 - alpha));
        end
    else
        if observationWeights
            shrinkFactor = weights*(X0.^2) + lam*(1 - alpha);
        else
            shrinkFactor = (1/N) * ones(1,N)*(X0.^2) + lam*(1 - alpha);
        end
    end

    % Iterative coordinate descent until converged
    for numIter = 1:maxIter
        
        bold = b;

        [b,active] = cdescentCycle(X0,wX0,Y0, ...
            b,active,totalweight,shrinkFactor,threshold);
        
        if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) < reltol
            % Cycling over the active set converged.
            % Do one full pass through the predictors.
            % If there is no predictor added to the active set, we're done.
            % Otherwise, resume the coordinate descent iterations.
            bold = b;
            potentially_active = thresholdScreen(X0,wX0,Y0,b,ever_active,threshold);
            if any(potentially_active)
                new_active = active | potentially_active;
                [b,new_active] = cdescentCycle(X0,wX0,Y0, ...
                    b,new_active,totalweight,shrinkFactor,threshold);
            else
                new_active = active;
            end

            if isequal(new_active, active)
                break
            else
                active = new_active;
            end
            
            if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) < reltol
                break
            end
        end

        if numIter == maxIter
            warning(message('stats:lasso:MaxIterReached',num2str(lam)));
        end
    end
    
    B(:,i) = b;
    
    % Halt if maximum model size ('DFmax') has been met or exceeded.
    if sum(active) > dfmax
        % truncate B and lambda output arguments
        lambda = lambda(1:(i-1));
        B = B(:,1:(i-1));
        break
    end
    
    % Halt if we have exceeded a threshold on the percent of
    % residual variance left unexplained.
    if ~userSuppliedLambda
        % Calculate mse of the current fit
        bsig = b ./ sigmaX';
        fit = [Intercept_Flag*ones(size(X,1),1) X] * [(muY-muX*bsig); bsig];
        residuals = bsxfun(@minus, Y, fit);
        if ~observationWeights
            mspe = mean(residuals.^2);
        else
            % This line relies on the weights having been normalized.
            mspe = weights * (residuals.^2);
        end
        if mspe < 1.0e-3 * nullMSE
            lambda = lambda(1:i);
            B = B(:,1:i);
            break
        end
    end
    
end % of lambda sequence

% ------------------------------------------
% Unwind the centering and scaling (if any)
% ------------------------------------------

B = bsxfun(@rdivide, B, sigmaX');
B(~ever_active,:) = 0;
Intercept = muY-muX*B;   % Intercept =TRUE

% ------------------------------------------
% Calculate Mean Prediction Squared Error
% ------------------------------------------

BwithI = [Intercept; B];
fits = [ones(size(X,1),1) X]*BwithI;
residuals = bsxfun(@minus, Y, fits);
if ~observationWeights
    mspe = mean(residuals.^2);
else
    % This line relies on the weights having been normalized.
    mspe = weights * (residuals.^2);
end

end %-lassoFit


% ===================================================
%                 cdescentCycle() 
% ===================================================

function [b,active] = cdescentCycle(X0, wX0, Y0, ...
    b, active, totalweight, shrinkFactor, threshold)
[b,active] = internal.stats.lassoCoordDescentCycle(X0,wX0,Y0,b,active',totalweight,shrinkFactor',threshold); 
active = active';
end %-cdescentCycle

% ===================================================
%                 thresholdScreen() 
% ===================================================

function potentially_active = thresholdScreen(X0, wX0, Y0, ...
    b, active, threshold)
r = internal.stats.lassoYminusXb( X0, Y0, b, active' );
potentially_active = abs(r' *wX0) > threshold;
end %-thresholdScreen

% ===================================================
%                 computeLambdaMaX() 
% ===================================================

function [lambdaMax, nullMSE] = computeLambdaMax(X, Y, weights, alpha, standardize)
%
% lambdaMax is the penalty term (lambda) beyond which coefficients
% are guaranteed to be all zero.
%
% nullMse is the mse of the fit using just a constant term.
% It is provided in this function as a convenience, because it needs 
% to be calculated in the same context as lambdaMax whenever
% lambdaMax is calculated.

if ~isempty(weights)
    observationWeights = true;
    weights = weights(:)';        
    % Normalized weights are used for standardization and calculating lambdaMax.
    normalizedweights = weights / sum(weights);
else
    observationWeights = false;
end

[N,~] = size(X);

% If we were asked to standardize the predictors, do so here because
% the calculation of lambdaMax needs the predictors as we will use
% them to perform fits.

if standardize
    % If X has any constant columns, we want to protect against
    % divide-by-zero in normalizing variances.
    constantPredictors = (range(X)==0);

    if ~observationWeights
        % Center and scale
        [X0,~,~] = zscore(X,1);
    else
        % Weighted center and scale
        muX = normalizedweights * X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = sqrt( normalizedweights * (X0.^2) );
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
        X0 = bsxfun(@rdivide, X0, sigmaX);
    end
else
    if ~observationWeights
        % Center
        muX = mean(X,1);
        X0 = bsxfun(@minus,X,muX);
    else
        % Weighted center
        muX = normalizedweights(:)' * X;
        X0 = bsxfun(@minus,X,muX);
    end
end

% If using observation weights, make a weighted copy of the 
% predictor matrix, for use in weighted dot products.

if observationWeights
    wX0 = bsxfun(@times, X0, weights');
end

if ~observationWeights
    muY = mean(Y);
else
    muY = weights*Y;
end
% Y0 = bsxfun(@minus,Y,muY);
Y0 = Y - muY;

% Calculate max lambda that permits non-zero coefficients
%
if ~observationWeights
    dotp = abs(X0' * Y0);
    lambdaMax = max(dotp) / (N*alpha);
else
    dotp = abs(sum(bsxfun(@times, wX0, Y0)));
    lambdaMax = max(dotp) / alpha;
end

if ~observationWeights
    nullMSE = mean(Y0.^2);
else
    % This works because weights are normalized and Y0 is already
    % weight-centered.
    nullMSE = weights * (Y0.^2);
end
end

