%% Generate n samples from a d-variate Gaussian distribution
%
% function X=mvgauss(mu, Sigma, n)
%
% Parameters:
%   mu    = mean [d x 1]
%   Sigma = variance-covariance matrix [d x d]
%
% Returns:
%   X     = n samples form a d-variate gaussian [n x d]
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2010
function X=mvgauss(mu, Sigma, n)

% d-variate normal
d=length(mu);

if(d == 1)      % if univariate normal requested
    X=normrnd(mu, sqrt(Sigma), n, 1);
   
else
    A=chol((Sigma+Sigma')./2);   % Cholesky decomposition, Sigma = A'A
    
    X=A'*randn(d,n);            % N(0,Sigma)
    if((mu'*mu)>0)              % if non zero mean
        X=bsxfun(@plus,X,mu);   % add the mean
    end;
    X=X';   % make sure X is [dxn]
end;

return;