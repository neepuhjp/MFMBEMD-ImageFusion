function [Q] = getTotalEvaluationValues(I_input1,I_input2,fuseImg)

bit1  = size(I_input1,3);
if bit1 > 1
    img1 = rgb2gray(I_input1);
else
     img1 = I_input1; 
end

bit2  = size(I_input2,3);
if bit2 > 1
    img2 = rgb2gray(I_input2);
else
     img2 = I_input2; 
end

bit3 = size(fuseImg,3);
if bit3 > 1
     imgf = uint8(rgb2gray(fuseImg));
else
     imgf = uint8(fuseImg); 
end


Q(1)=metricMI(img1,img2,imgf,1);
Q(2)=analysis_fmi(img1,img2,imgf); 
Q(3)=metricWang(img1,img2,imgf);
Q(4)=metricXydeas(img1,img2,imgf); %Q_abf
Q(5)=metricZhao(img1,img2,imgf);
Q(6)=metricPeilla(img1,img2,imgf,3);  %
Q(7)=metricYang(img1,img2,imgf); %Q_xyf
Q(8)=metricChenBlum(img1,img2,imgf);
Q(9)=metricChen(img1,img2,imgf);   %
Q(10)=vifp_mscale(img1,imgf) + vifp_mscale(img2,imgf);
end

