function SF = SF_evaluation(MF)
RF = diff(MF,1,1);
RF = sqrt(mean(mean(RF.^2)));
CF = diff(MF,1,2);
CF = sqrt(mean(mean(CF.^2)));
SF = sqrt(RF^2+CF^2);
end
