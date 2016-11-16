function [rateMat]=rate(matIn,timeRes)
% Calculates rate of apical constriction and myosin accumulation for all 
% cells in a data set.
% 
% It uses the avg of t+1 and t-1 rates.
% 
% **For future, may want to invert constriction rates such that correlation 
% is positive...**
% 
% INPUT: smoothed data(dataArea or dataMyo)
%        timeRes (seconds per stack, during image acquisition)
% OUTPUT: matrix that contain appropriate rate

rateMat=zeros(size(matIn));

for cell_index=1:size(matIn,2)
    for n=1:size(matIn,1)
        if (n==1)
            rateMat(n,cell_index)=(matIn(n+1,cell_index)-matIn(n,cell_index))/timeRes;        
        elseif (n==size(matIn,1))
            rateMat(n,cell_index)=(matIn(n,cell_index)-matIn(n-1,cell_index))/timeRes;
        else
            rateMat(n,cell_index)= (matIn(n+1,cell_index)-matIn(n-1,cell_index))/(2*timeRes);
        end      
    end
end

