 
function CI=CI_Sun_Normal_PriorPPI(Summary,sig_level)
alpha=sig_level/2;
% x = data;                      % Create Data
SEM = Summary.Std;  %/sqrt(length(x));  % Standard Error  % std has been normanized by sqrt(N-1).    
ts = tinv([alpha  1-alpha],1e4-1);      % T-Score
CI = Summary.Mean + ts.*SEM;                      % Confidence Intervals
CI = CI(1:size(CI,1),:);

