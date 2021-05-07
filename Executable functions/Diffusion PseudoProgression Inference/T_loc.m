function [T1, phi0] = T_loc(data, nsig, W)
%performs locally scaled diffusion map, suitable for rather small data
%size such as clinical transcriptomic data


% data=high dimensional matrix data to be analysied 
% nsig=kernel width will be locally adjusted to the distance to the nsig'th nearest neighbor
% W=weight coefficients by clinical information
% T1 : transition matrix of data
% phi0: zeroth eig.vector of T1, corresponding to eig.val=1


tic

n1=size(data,1);
d2=pdist(data,'euclidean').^2;
d2=squareform(d2);

[idnn_s,dists]=knnsearch(data,data,nsig);
% [idnn_s,dists]=knnsearch(data);
sigma=dists(:,nsig);

% S=sigma*sigma';
S2 = bsxfun(@plus,sigma.^2,sigma.^2');  
Sw=S2./W;
% W1=sqrt(2*S./S2).*exp(-dw./S2);
W1=exp(-d2./Sw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D1=sum(W1,2);

q=bsxfun(@times,reshape(D1,1,n1),reshape(D1,n1,1)).^1;

W1=W1./q;
%W1(1:n1+1:end)=0; % not necessary with local sigma
W1(d2==0)=0; % not necessary with local sigma

D1_=diag(sum(W1,2));
sum(sum(isnan(D1)))
sum(sum(isnan(W1)))
T1=D1_^(-0.5)*W1*D1_^(-0.5);
phi0 = diag(D1_)/sqrt(sum(diag(D1_).^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%.04f', toc/60);
    

    
