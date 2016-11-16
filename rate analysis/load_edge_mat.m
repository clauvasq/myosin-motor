function [embryoOut]=load_edge_mat(embryoInfo)
% This function takes the input embryoInfo, a cell array with the following
% inputs {name,type,timeRes,t-initial,t-final} 
% where t-initial and t-final are times that are segmented.
% it outputs a structure containing the input info along with myosin, area
% and anisotropy data.

for i=1:size(embryoInfo,1)
   dataOut(i).name=embryoInfo{i,1};
   dataOut(i).type=embryoInfo{i,2};
   dataOut(i).timeRes=embryoInfo{i,3};
   
   % where measurements folders live, concat folder name of each specific
   % embryo Dropbox' (MIT)'/
   startFilename=strcat('~/Documents/MATLAB/EDGE-output/', embryoInfo{i,1});
   filenames={strcat(startFilename,'/Myosin--myosin_intensity--Myosin intensity.mat'),...
       strcat(startFilename,'/Membranes--basic_2d--Area.mat'),...
       strcat(startFilename,'/Membranes--ellipse_properties--Anisotropy-xy.mat'),...
       strcat(startFilename,'/Membranes--ellipse_properties--Anisotropy.mat')};
   
    z = 1;
   load(filenames{1}) 
   dataOut(i).myo=squeeze(cell2mat(data(:,z,:)));
   clear data name unit
   
   load(filenames{2}) 
   dataOut(i).area=squeeze(cell2mat(data(:,z,:)));
   clear data name unit
   
   load(filenames{3}) 
   dataOut(i).anisotropyXY=squeeze(cell2mat(data(:,z,:)));
   clear data name unit
   
   load(filenames{4}) 
   dataOut(i).anisotropy=squeeze(cell2mat(data(:,z,:)));
   clear data name unit
   
   if size(dataOut(i).myo,1)~=size(embryoInfo{i,4}:embryoInfo{i,5},2)
        dataOut(i).myo=dataOut(i).myo(embryoInfo{i,4}:embryoInfo{i,5},:);
        dataOut(i).area=dataOut(i).area(embryoInfo{i,4}:embryoInfo{i,5},:);
        dataOut(i).anisotropyXY=dataOut(i).anisotropyXY(embryoInfo{i,4}:embryoInfo{i,5},:);
        dataOut(i).anisotropy=dataOut(i).anisotropy(embryoInfo{i,4}:embryoInfo{i,5},:);
   end
   dataOut(i).timeSec=dataOut(i).timeRes*[0:size(dataOut(i).myo,1)-1];
   dataOut(i).timeMin=dataOut(i).timeSec./60;
end
embryoOut=dataOut;
end