%data=genes;

function [rankeddata]=rankdata(data)
tic
rankeddata=zeros(size(data,1),size(data,2));

for g=1:size(data,2)
    A=unique(sort(data(:,g)));
    for sample=1:size(data,1)
       
        rankeddata(sample,g)=find(data(sample,g)==A);
    end
end
  
toc/60
    