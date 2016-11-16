%% Uses kymogrpahs made in Fiji using kymograph_max_thr_binary macro
% in short makes kymograph that is 100pix20pi along center of cut
% makes maximum projection of 20 kymographs
% thresholds and then makes image binary
% this code will import that data, find the distance between the tissue
% displacements on each side of the cut, and plots that

%% Load filenames
input_pw = '~/Documents/Martin Lab/Projects/sqhXX_2P/sqhXX_2P vline/kymograph_analysis/kymo_binary_using/';
kymoFiles = dir(strcat(input_pw));
%removes hidden filename
kymoFiles = kymoFiles(arrayfun(@(x) ~strcmp(x.name(1),'.'),kymoFiles));
numFiles = numel(kymoFiles);
res = 0.2; % spatial res. um/pix
timeRes = 1.3; %temporal res sec/stack
tAblation = 10; %frame before ablation
%% Open each file and do stuff to them
% rows = y-dim of kymograph (originally x-dimension of cutting image)
% columns = time-dimension
% zero = black, one = white
cutLoc = 50;
for i = 1:numFiles
    filename = strcat(input_pw, kymoFiles(i).name);
    rawData = double(imread(filename));
    rawData = rawData./255;  
%     clean data to remove outlier points using bwmorph
    cleanData = bwmorph(rawData,'clean');
    cleanData = bwmorph(rawData,'fill');
    cleanData = bwmorph(rawData, 'majority');    
%     find consequtive runs of ones
    cutIndex = NaN(1,size(cleanData,2));
    cutSpan = NaN(1,size(cleanData,2));
    for j = 1:size(cleanData,2)
        a = cleanData(:,j);
        dsig = diff([0 a.' 0]);
        startIndices = find(dsig > 0); 
        endIndices = find(dsig < 0)-1;
        b = sum(startIndices <= cutLoc);
        if b > 0
            cutIndex(j) = startIndices(b);
            cutSpan(j) = endIndices(b) - startIndices(b)+1;
        else
            cutIndex(j) = cutLoc;
            cutSpan(j) = 2; 
        end
    end
%     save appropriate variables in a structure
    emb(i).rawMat = rawData;
    emb(i).cleanMat = cleanData;
    emb(i).cutIndex = cutIndex;
    emb(i).cutSpan = cutSpan;
    emb(i).cutSpanUm = cutSpan*res;
end

%% Plotting (plot each cut as an individual line)
% i = 1-10 is AE; 11-21 is TA; 22-36 is TS
embAE = [1, 10];
embTA = [11, 21];
embTS = [22, 36];

figure
for i = 1:3
    if i == 1
        startEmb = embTS(1);
        endEmb = embTS(2);
    elseif i == 2
        startEmb = embTA(1);
        endEmb = embTA(2);
    else
        startEmb = embAE(1);
        endEmb = embAE(2);
    end
    subplot(1,3,i)
    hold on
    for j = startEmb:endEmb
        emb(j).timeMat = ((1:size(emb(j).cutSpan,2)) - tAblation).*timeRes; 
        plot(emb(j).timeMat, emb(j).cutSpanUm)
    end
    hold off
end

%% make matrices that have rows for each embryo with only frames 1 - 50

endFrame = 50;
for i = 1:3
    if i == 1
        startEmb = embTS(1);
        endEmb = embTS(2);        
    elseif i == 2
        startEmb = embTA(1);
        endEmb = embTA(2);
    else
        startEmb = embAE(1);
        endEmb = embAE(2);
    end
    numEmb = endEmb - startEmb + 1;
    dispTemp = NaN(numEmb,endFrame);
    for j = startEmb:endEmb
        index_j = j - startEmb + 1;
        dispTemp(index_j, :) = emb(j).cutSpanUm(1:endFrame); 
    end
%     make such that disp = 0 at tAblation
    dispTemp = bsxfun(@minus, dispTemp, dispTemp(:,tAblation));
    if i == 1
        dispTS = dispTemp./2;
    elseif i == 2
        dispTA = dispTemp./2;
    else
        dispAE = dispTemp./2;
    end
end

%% Plot averages
timeMat = ([1:50] - tAblation).*timeRes;
figure
% plot(timeMat, mean(dispTS), timeMat, mean(dispTA), timeMat, mean(dispAE))
hold on
shadedErrorBar(timeMat,dispTS,{@mean,@std},'c')
shadedErrorBar(timeMat,dispTA,{@mean,@std},'b')
shadedErrorBar(timeMat,dispAE,{@mean,@std},'m')
hold off
