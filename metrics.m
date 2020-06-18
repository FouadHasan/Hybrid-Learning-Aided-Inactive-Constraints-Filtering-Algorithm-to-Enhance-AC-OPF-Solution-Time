clc
clear all


x=597;%put no of bus(for voltage constraints) or line(for line constraints)

TN=1752000;
FN=0;
FP=26335;
TP=12665;

%no of scenario
s=(TP+TN+FP+FN)/x
%%
c=s*x*2;
%%
TN=c-(TP+FP+FN);

Accuracy=(TP+TN)/(TP+TN+FP+FN)
Misclassification=(FP+FN)/(TP+TN+FP+FN)
FNP=(FN)/(TP+TN+FP+FN);
FPP=(FP)/(TP+TN+FP+FN);
TNP=(TN)/(TP+TN+FP+FN);
TPP=(TP)/(TP+TN+FP+FN);
disp('-------------------')

NPV=TN/(FN+TN);
% FOR=1-NPV 
PPV=TP/(FP+TP);          
% FDR=1-PPV                    
TPR=(TP)/(TP+FN);
% FNR=1-TPR                   
TNR=TN/(TN+FP);
% FPR=1-TNR                      

ALL_Metrics=[FNP,FPP,TNP,TPP,NPV,PPV,TPR,TNR,Misclassification,Accuracy]*100       

