function [fuseImg] = SingleImgPatchFusion(I_input1, I_input2, decNum, winSize, stepLen,imageType)

[mrows, ncols,bit]  = size(I_input1);
fuseImg = I_input1;

for i = 1:bit
    %Construct a multichannel signal from the input images to be fused
    InImgs = zeros(mrows,ncols,2);
    InImgs(:,:,1) = I_input1(:,:,i);
    InImgs(:,:,2) = I_input2(:,:,i);
    
    [decRlt,~] = multiChannelMFEMD(InImgs,decNum);
    [IMF1,residue1,IMF2,residue2] = prepareMFMEMDFusionInput(decRlt,mrows, ncols);
    
    [fuseIMFs,fuseRes] = patchFusionViaEMD(IMF1,residue1,IMF2,residue2,winSize, stepLen,imageType);
    
    fuseImgIMFs(:,:,:,i) = fuseIMFs;
    fuseImgRes(:,:,i) = fuseRes;
end

%reconstruction
for i = 1:bit
    fuseImg(:,:,i) = sum(fuseImgIMFs(:,:,:,i),3) + double(fuseImgRes(:,:,i));
end

end



function [IMF1,residue1,IMF2,residue2] = prepareMFMEMDFusionInput(decRlt,mrows, ncols)

IMFNum = size(decRlt,2)-1;
IMF1 = zeros(mrows,ncols,IMFNum);
IMF2 = zeros(mrows,ncols,IMFNum);
for m = 1:IMFNum
    tempCIMFs = decRlt{1,m};
    IMF1(:,:,m) = tempCIMFs(:,:,1);
    IMF2(:,:,m) = tempCIMFs(:,:,2);
end

residue = decRlt{1,end};
residue1 = residue(:,:,1);     residue2 = residue(:,:,2);

end


function [fuseIMFs,fuseRes] = patchFusionViaEMD(IMFs1,res1,IMFs2,res2,winSize,stepLen,imageType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = size(IMFs1,1);
h = size(IMFs1,2);
IMFNum = size(IMFs1,3);
fuseIMFs = zeros(w,h,IMFNum);
fuseRes = zeros(w,h);

countFun = zeros(w, h);
countWindow = ones(winSize, winSize);

xFilterMax = w-winSize+1;
yFilterMax = h-winSize+1;

[IMFWeights,resWeights] = computeFusionWeightsViaEMD(IMFs1,IMFs2,winSize,imageType);

sampX = 1 : stepLen : xFilterMax;
sampX = [sampX xFilterMax];
sampY = 1 : stepLen : yFilterMax;
sampY = [sampY  yFilterMax];
offset = winSize-1;

for row = 1 : length(sampX)
    for col = 1 : length(sampY)
        i = sampX(row);
        j = sampY(col);
        blocks1 = IMFs1(i:i+offset, j:j+offset, :); 
        blocks2 = IMFs2(i:i+offset, j:j+offset, :);
        IMFBlock = zeros(winSize,winSize,IMFNum);
        
        for index = 1:IMFNum
           IMFBlock(:,:,index) = IMFWeights(i,j,index,1) * blocks1(:,:,index) + IMFWeights(i,j,index,2)  * blocks2(:,:,index);
        end
        fuseIMFs(i:i+offset, j:j+offset, :) = fuseIMFs(i:i+offset, j:j+offset, :) + IMFBlock;
          
        blocks1 = res1(i:i+offset, j:j+offset, :); 
        blocks2 = res2(i:i+offset, j:j+offset, :); 
        
        if imageType == -1 % the fusion for the residue of the multimodal images
            blocks_mean1 = mean(mean(blocks1));
            blocks_mean2 = mean(mean(blocks2));
            p = 6;
            y1 = blocks_mean1^p;
            y2 = blocks_mean2^p;
            ResBlock = resWeights(i,j,1)  * (blocks1 - blocks_mean1) + resWeights(i,j,2) * (blocks2 - blocks_mean2) + y1 / (y1 + y2) * blocks_mean1 + y2 / (y1 + y2) *  blocks_mean2;
        else   % the fusion for the residue of the multi-focus images
            ResBlock = resWeights(i,j,1)  * blocks1 + resWeights(i,j,2) * blocks2;
        end
        
        fuseRes(i:i+offset, j:j+offset, :) = fuseRes(i:i+offset, j:j+offset, :) + ResBlock;                 
        countFun(i:i+offset, j:j+offset, :) = countFun(i:i+offset, j:j+offset, :) + countWindow;
    end
end

fuseIMFs = fuseIMFs ./ countFun;
fuseRes = fuseRes ./ countFun;
end

function [IMFWeights,resWeights] = computeFusionWeightsViaEMD(IMFs1,IMFs2,winSize,imageType)

w = size(IMFs1,1);
h = size(IMFs1,2);
IMFNum = size(IMFs1,3);
xFilterMax = w-winSize+1;
yFilterMax = h-winSize+1;

window = ones(winSize);
window = window / sum(window(:));

IMFEnergy = zeros(xFilterMax, yFilterMax, IMFNum, 2);
IMFWeights = zeros(xFilterMax, yFilterMax, IMFNum, 2);

for i = 1 : IMFNum
    tempImg = filter2(window, IMFs1(:, :, i) .* IMFs1(:, :, i), 'valid');
    IMFEnergy(:,:,i,1) = tempImg;
    tempImg = filter2(window, IMFs2(:, :, i) .* IMFs2(:, :, i), 'valid');
    IMFEnergy(:,:,i,2) = tempImg;
end

[~, Labels] = max(IMFEnergy,[],4);

for i = 1:2
    tempIndex = zeros(xFilterMax,yFilterMax,IMFNum);
    tempIndex(Labels == i) = 1;
    IMFWeights(:,:,:,i) = tempIndex;
end

if imageType == 1  % multi-focus images
    resWeights(:,:,1) = IMFWeights(:,:,1,1);
    resWeights(:,:,2) = IMFWeights(:,:,1,2);
else   %% multi-modal images   
    p = 6;%
    resWeights =  zeros(xFilterMax, yFilterMax, 2);
    TotalIMFEnergy = zeros(xFilterMax, yFilterMax, 2);
    for i = 1 : IMFNum
        TotalIMFEnergy(:,:,1) = TotalIMFEnergy(:,:,1) + IMFEnergy(:,:,i,1);
        TotalIMFEnergy(:,:,2) = TotalIMFEnergy(:,:,2) + IMFEnergy(:,:,i,2);
    end
    
    TotalIMFEnergy = TotalIMFEnergy.^p;
    resWeights(:,:,1) = TotalIMFEnergy(:,:,1) ./ (TotalIMFEnergy(:,:,1) + TotalIMFEnergy(:,:,2)+ 0.000001);
    resWeights(:,:,2) = TotalIMFEnergy(:,:,2) ./ (TotalIMFEnergy(:,:,1) + TotalIMFEnergy(:,:,2)+ 0.000001);
end

end









