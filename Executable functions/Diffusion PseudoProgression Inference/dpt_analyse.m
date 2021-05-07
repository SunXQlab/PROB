function DPT = dpt_analyse(M)
% Given the accumulated transition matrix M, dpt_analyse performs
% pseudotime ordering 
% M: the accumulated transition matrix

% DPT is the diffusion pseudotime in respect to the root cell





    n=size(M,1);

%   
%      dptbi=zeros(1,n);
     dptbi=dpt_to_root(M,1);

     dpt_branch{1}=mean(dptbi,1);

     DPT=dpt_branch{1}';

end

