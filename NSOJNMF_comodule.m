function [Co_module, Subpattern1, Subpattern2, Subpattern3] = NSOJNMF_comodule(X1,X2,X3,W,H1,H2,H3,tt0,tt1,tt2,tt3) 
% 
% Compute the mean and standard deviation of each rows in H to determine the module member. 
%
% INPUT:
% X1, X2, X3 : raw matrices
% W, H1, H2, H3 : obtained from NSOJNMF
% tt*: thresholds
%
% OUTPUT:
%
% Co_module is the list of sample and mRNAs, lncRNAs, miRNAs in a Co-module.
% Subpattern1, Subpattern2 and  Subpattern3 is the corresponding sub-matrices for each Co-module.

m1 = size(H1,2);
m2 = size(H2,2);
m3 = size(H3,2);
n = size(W,1);
K = size(W,2);

MW = mean(W,1);
MH1 = mean(H1,2); 
MH2 = mean(H2,2); 
MH3 = mean(H3,2);

VW = std(W,0,1); % std of each col
VH1 = std(H1,0,2); % std of each row
VH2 = std(H2,0,2);
VH3 = std(H3,0,2);

% Co-Module
for i = 1:K
    %W
    r = find(W(:,i)> MW(i) + tt0*VW(i));
    Co_module{i,1} = r';
    
    %H1
    c1 = find(H1(i,:)> MH1(i) + tt1*VH1(i));
    Co_module{i,2} = c1'; 
    Subpattern1{i} = X1(r,c1');
    
    %H2
    c2 = find(H2(i,:)> MH2(i) + tt2*VH2(i));
    Co_module{i,3} = c2'; 
    Subpattern2{i} = X2(r,c2');
    
    %H3
    c3 = find(H3(i,:)> MH3(i) + tt3*VH3(i));
    Co_module{i,4} = c3'; 
    Subpattern3{i} = X3(r,c3');
  
end
        
end

