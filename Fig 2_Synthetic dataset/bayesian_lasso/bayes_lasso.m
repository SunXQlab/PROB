%% Implementation of the Bayesian LASSO linear regression method.
%  Algorithm by Park and Casella, 2008.
%
% function [beta, sigma2, tau2, lambda]=bayes_lasso(X,y)
%
% Parameters:
%   X     = centered regressor matrix [n x p]
%   y     = centered response vector [n x 1]
%   niter = number of iterations for the Gibbs sampler [optional,default=1e5]
%
% Returns:
%   beta   = shrunken regression parameters [p x niter]
%   sigma2 = noise variance  [1 x niter]
%   tau2   = prior variance for the regressors [p x niter]
%   lambda = exponential hyperparameter [1 x niter]
%
% References:
%
% The Bayesian Lasso
% Trevor Park and George Casella
% Journal of the American Statistical Association, Vol. 103, No. 482, pp.
% 681-686, 2008.
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2008
function [beta, sigma2, tau2, lambda]=bayes_lasso(X,y,niter)

% default niter
if(nargin < 3)
    niter=1e5;
end;

%% Initialize parameters
[n,k]=size(X);              % size of the data set

XtX=X'*X;                   % precompute for speed
Xy=X'*y;        

% store Gibbs trace for each parameter
beta=zeros(k,niter);      
sigma2=zeros(1,niter);
tau2=zeros(k,niter);
lambda=zeros(1,niter);

%% Initial Gibbs step
iter=1;
beta(:,iter)=X\y;                                           % set beta to the LS estimates

e=y-X*beta(:,iter);
sigma2(1,iter)=e'*e/n;                                      % set sigma2 to the LS estimate
tau2(:,iter)=beta(:,iter).^2;                               % set tau2 to the LS estimates (using beta as 'data')
lambda(1,iter)=k*sqrt( sigma2(1,iter) ) / sum(abs( beta(:,iter) ) );    % initialize the hyperparameter

%% Gibbs sampler
iter=2;
while(iter <= niter)
    
    %% debug info
    if(mod(iter,100) == 0)
        fprintf('Iter(%d)\n', iter);
    end;
    
    %% Sample from the beta conditional dist.
    Dinv=diag(1./tau2(:,iter-1));
    A=(XtX + Dinv);
    mu=A\Xy;                     
    Sigma=sigma2(1,iter-1) * (A \ eye(size(A)));
    beta(:,iter)=mvgauss(mu, Sigma, 1)';   % beta | (.) ~ N(mu, Sigma)
    
    %% Sample sigma2 from the conditional dist.
    shape=(n-1)/2 + k/2;
    e=y-X*beta(:,iter);
    scale=e'*e/2 + beta(:,iter)'*Dinv*beta(:,iter)/2;
    sigma2(1,iter)=1/gamrnd(shape, 1/scale);
    
    %% Sample tau2 from the conditional
    mu_hat=sqrt( (lambda(1,iter-1)^2*sigma2(1,iter))./ beta(:,iter).^2 );
    lambda_hat=lambda(1,iter-1)^2;
    for i=1:k,
        tau2(i,iter)=1/randraw('ig', [mu_hat(i),lambda_hat]);
    end;
    
    %% Empirical Bayes Update lambda
    if mod(iter,10)==0,
        low=iter-9;
        high=iter;
        lambda(1,iter)=sqrt( (2*k) / sum(mean(tau2(:,low:high),2)));
    else
        lambda(1,iter)=lambda(1,iter-1);
    end;        
    
    iter=iter+1;    % next iteration
end;

return;
