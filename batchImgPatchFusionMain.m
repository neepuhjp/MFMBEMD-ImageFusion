clear;
% close all;

addpath(genpath(pwd)); 


% inputfile = '.\inputImages\grayscale\*.bmp'; 
% inputfile = '.\inputImages\color\*.bmp'; 
inputfile = '.\inputImages\medical\*.bmp'; 
%  inputfile = '.\inputImages\visible-infrared\*.bmp'; 

namelist = dir(inputfile);
outputfile = '.\outputImages\';

decNum = 2;   %Decomposition number of MBEMD, multifocus images: 1; multimodal images: 2  
winSize = 15; %Window size, gray-scale mutlifocus: 31; color multifocus: 15; medical multimodal: 15; infrared-visible:33
stepLen = round(winSize - (winSize - 2)); %the parameter related to the overlapping size, multifocus: round(winSize - 1.0/6 * winSize); multimodal: (winSize - (winSize - 2))£»
imageType = -1; %Image type, mutlifocus: 1; multimodal: -1

len = length(namelist);
TotalEvalResult = zeros(len/2,10);

T = {'Q_MI','Q_FMI','Q_NCIE','Q_G','Q_P','Q_E','Q_Y','Q_CB','Q_CV','Q_VIF'};
B1 = {'inputName'}; B2 = {'outputName'};

for i = 1:len/2
   inputName1 = namelist(2*i-1).name; inputName2 = namelist(2*i).name;
   I_input1 = imread(inputName1);
   I_input2 = imread(inputName2);   
   
   tic; 
   [output] = SingleImgPatchFusion(I_input1, I_input2, decNum,winSize,stepLen,imageType); 
   toc;
   
   methodName = 'MFMEMD';
   outputName = [outputfile,strcat(methodName,[inputName1(1:end-4),'.bmp'])];
   imwrite(uint8(output),outputName);    
   
   outputName = strcat(methodName,[inputName1(1:end-4),'.bmp']);
   B1 = [B1;inputName1(1:end-4)];
   B2 = [B2;outputName(1:end-4)];
   TotalEvalResult(i,:) = getTotalEvaluationValues(I_input1,I_input2,uint8(output));
end

A = [T; num2cell(TotalEvalResult)];

xlswrite(strcat(methodName,'TotalEvalResult.xlsx'),cell(100,100)); 
xlswrite(strcat(methodName,'TotalEvalResult.xlsx'),[B1,B2,A]);

EvalResult(1,:) = mean(TotalEvalResult,1); 
EvalResult(2,:) = std(TotalEvalResult,0,1);
C = [T;num2cell(EvalResult)];
B = {'Name';'Mean';'std'};
xlswrite(strcat(methodName,'meanStdEvalResult.xlsx'),cell(100,100)); 
xlswrite(strcat(methodName,'meanStdEvalResult.xlsx'),[B,C]);



 
