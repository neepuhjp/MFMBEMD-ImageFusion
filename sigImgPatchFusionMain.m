
addpath(genpath(pwd)); 

clear;
[f, pathname,FilterIndex] = uigetfile({'*.*','*.allFiles'},'Open an image');
file = [pathname,f];

if FilterIndex == 1
   I_input1 = imread(file);
else
    error('There is a wrong in opening File!');
end

[row1, col1,bit1] = size(I_input1);



[f, pathname,FilterIndex] = uigetfile({'*.*','*.allFiles'},'Open an image');
file = [pathname,f];

if FilterIndex == 1
   I_input2 = imread(file);
else
    error('There is a wrong in opening File!');
end


[row2, col2,bit2] = size(I_input2);

decNum = 2;   %the decomposition number of MBEMD, multi-focus images: 1; multi-modal images: 2  
winSize = 33; %the window size, gray mutli-focus: 31; color multi-focus: 15; medical multi-modal: 15; infrared-visible:33
stepLen = round(winSize - (winSize - 2)); %the parameter related to the overlapping size, multifocus: round(winSize - 1.0/6 * winSize); multimodal: (winSize - (winSize - 2))£»
imageType = -1; %image type, mutli-focus: 1; multimodal: -1


tic;
[output] = SingleImgPatchFusion(I_input1, I_input2, decNum,winSize,stepLen,imageType);  
toc;


methodName = 'MFMEMD';
outputName = strcat(methodName,[f(1:end-4),'.bmp']);
imwrite(uint8(output),outputName);

EvalResult = getTotalEvaluationValues(I_input1,I_input2,uint8(output));

figure; imshow(I_input1);
figure; imshow(I_input2);
figure;imshow(uint8(output)); 






