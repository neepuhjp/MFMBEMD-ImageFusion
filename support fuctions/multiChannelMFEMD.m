function [decRlt,avgDis] = multiChannelMFEMD(mchImg,numIMF,avgDisTh,extrNumTh)
%This function implements the decomposition of a multi-channel image (e.g., rgb images) by the multi-channel bidimensional EMD based on
%morphological filter, which is a multi-channel extension of the enhanced ast empirical mode decomposition method
%(M. Trusiak, M. Wielgus, K. Patorski, Advanced processing of optical fringe patterns by automated selective reconstruction and enhanced fast empirical mode decomposition, Optics & Lasers in Engineering 52 (2014) 230{240.)

%mchImg - the input image to be decomposed, which can be a rgb image
%numIMF - number of decomposition level
%avgDisTh - the threshold of the mean distance between lenrema, and the default value is half of the maximum of the image height and height.
%extrNumTh: the minimum of extremum number accepted in the decomposition.
%decRlt£º the decomposition results including IMFs and the residue.
%avgDis: the mean extremum distance of each residue during the decomposition

% example: decResult = multiChannelMFEMD(I_input,10,max(row,col)/2,3);


[height,width,dim] = size(mchImg);

if nargin < 3
    avgDisTh = max(height,width) / 2.0;
    extrNumTh = 0;
elseif nargin < 4
    extrNumTh = 0;
end


%-extend the image boundary
len = round(min(height,width)/3); %
newMchImg = zeros(height + 2*len, width + 2*len,dim);
for i = 1: dim
    Im = mchImg(:,:,i);
    p1 = flipud(Im(1:len,:)) ;
    p2 = flipud(Im(height-len+1:height,:));
    p3 = fliplr(Im(:,1:len));
    p4 = fliplr(Im(:,width-len+1:width)) ;
    p5 = Im(1:len,1:len);
    p6 = Im(1:len,width-len+1:width);
    p7 = Im(height-len+1:height,1:len);
    p8 = Im(height-len+1:height,width-len+1:width);
    Im = [fliplr(rot90(p5)), p1, fliplr(rot90(p6')); p3, Im, p4; fliplr(rot90(p7')), p2, fliplr(rot90(p8'))];
    newMchImg(:,:,i) = Im;
end

% imshow(uint8(newMchImg));

num = 1; % the number of the extracted IMF

[resMinAvDist, curExtrNum, ~] =  GetExAvDistViaMorph(newMchImg);

if curExtrNum <= extrNumTh
    decRlt{num} = newMchImg(len+1:height+len ,len+1:width+len,:);
    return;
end

minAvDistNew = resMinAvDist;
minAvDist = 0;
avgDis(num) = resMinAvDist;

while (num <= numIMF) %sifting processing
    %---Compute the size of morphological filters
    if minAvDistNew > minAvDist
        minAvDist = minAvDistNew;
    else
        minAvDist = round(2*minAvDist);
    end
    %---------------------------------------------------------------
    
    %---Stop the sifting procesing
    if resMinAvDist >= avgDisTh || curExtrNum <= extrNumTh
        decRlt{num} = newMchImg(len+1:height+len ,len+1:width+len,:);
        return;
    end
    
    [IMF,residue] = GetIMF(newMchImg,minAvDist);
    decRlt{num} = IMF(len+1:height+len ,len+1:width+len,:); %saving the residue
    newMchImg = residue;
    
    [resMinAvDist,curExtrNum, ~] =  GetExAvDistViaMorph(newMchImg);
    
    minAvDistNew = resMinAvDist;
    num = num + 1;
    avgDis(num) = resMinAvDist;
end

decRlt{num} = newMchImg(len+1:height+len ,len+1:width+len,:);
end

function [minAvDist,curExtrNum,avgDist] = GetExAvDistViaMorph(mchImg)
[M,N,dim] = size(mchImg);
filterMask = strel('square',3);
tempVec = [];
for i = 1:dim
    Img = mchImg(:,:,i);
    MaxEnv = imdilate(Img,filterMask);
    MinEnv = imerode(Img,filterMask);
    
    MaxMap = ~(Img - MaxEnv);
    MinMap = ~(Img - MinEnv);
    
    maxNum = sum(sum(MaxMap));
    minNum = sum(sum(MinMap));
    
    tempVec = [tempVec,maxNum,minNum] ;
    
    NumbExtrema = round(0.5*(maxNum + minNum));
    avgDist(i) = round(sqrt(N*M/NumbExtrema));
end

curExtrNum = min(tempVec);
minAvDist = min(avgDist);

if minAvDist == 1 %%%%%%%%%%%%%% the minimal average distance should be larger than 1
    minAvDist = 2;
end
end


function [IMF,residue] = GetIMF(mchImg,avgDist)

[M,N,dim] = size(mchImg);

maskDil = strel('square',avgDist);
maskConv = fspecial('average',[avgDist,avgDist]);

IMF = zeros(M,N,dim);
residue = zeros(M,N,dim);

for i = 1:dim
    Img = mchImg(:,:,i);
    MaxEnv = imdilate(Img,maskDil);
    MinEnv = -imdilate(-Img,maskDil);
    
    residue(:,:,i) = filter2(maskConv,0.5*(MaxEnv + MinEnv));
    IMF(:,:,i) = double(Img) - residue(:,:,i);
end

end

 






