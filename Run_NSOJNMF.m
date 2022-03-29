%
clear;clc;
cd('E:/NSOJNMF')

%%input data
mRNA = importdata('mRNA_log.txt');
lncRNA = importdata('lncRNA_log.txt');
miRNA = importdata('miRNA_log.txt');
theta1 = importdata('m_mRNAInte.txt');
theta2 = importdata('lnc_lncRNAInte.txt');
theta3 = importdata('mi_miRNAInte.txt');
R12 = importdata('m_lncRNAInte.txt');
R13 = importdata('mi_mRNAInte.txt');
R23 = importdata('mi_lncRNAInte.txt');

%extract samples and labes
mRNAs = mRNA.textdata;
lncRNAs = lncRNA.textdata;
miRNAs = miRNA.textdata;

%extract dataMatrix
X1 = mRNA.data;
X2 = lncRNA.data;
X3 = miRNA.data;
theta1 = theta1.data;
theta2 = theta2.data;
theta3 = theta3.data;
R12 = R12.data;
R13 = R13.data;
R23 = R23.data;
 
%% Run SOJNMF algorithm
[n,m1] = size(X1);
m2 = size(X2,2);
m3 = size(X3,2);
%parameters
nloop = 100; % 100; 
K = 200;
lambda1 = 0.001;
lambda2 = 0.001;
alpha = 0.01;
gamma = 10;
maxiter = 500; 
speak = 1;

bestW = zeros(n,K);
bestH1 = zeros(K,m1);
bestH2 = zeros(K,m2);
bestH3 = zeros(K,m3);

bestobj1=1000000000;
bestobj2=1000000000;
bestobj3=1000000000;

fid = fopen([ 'Record_K=' int2str(K) '_lambda1=' num2str(lambda1) '_lambda2=' num2str(lambda2) '_alpha=' num2str(alpha) 'gamma=' num2str(gamma) '.txt'],'wt+');

for iloop=1:nloop

    if speak 
        fprintf(1,' iteration %d\n',iloop); 
    end 
    
    [W,H1,H2,H3] = NSOJNMF(X1,X2,X3,theta1,theta2,theta3,R12,R13,R23,K,lambda1,lambda2,alpha,gamma,maxiter,speak,fid,iloop);

    
    % compute residue
    newobj1 = sum(sum((X1-W*H1).^2));
    newobj2 = sum(sum((X2-W*H2).^2));
    newobj3 = sum(sum((X3-W*H3).^2));
    
    if (newobj1<bestobj1)||(newobj2<bestobj2)||(newobj3<bestobj3)       
        bestobj1 = newobj1;
        bestobj2 = newobj2;
        bestobj3 = newobj3;
        bestW = W;
        bestH1 = H1;
        bestH2 = H2;
        bestH3 = H3;
    end
end
fclose(fid);

% compute the modules according to bestW, bestH1 and bestH2
W = bestW; 
H1 = bestH1; 
H2 = bestH2; 
H3 = bestH3;

% Output rule
% clear redundant data
clearvars -except X1 X2 X3 W H1 H2 H3 K isdouble mRNAs lncRNAs miRNAs
save NSOJNMF_2234Results.mat

tt0 = 3.5; tt1 = 3.5; tt2 = 3.5; tt3 = 3.5;
[Co_module, Subpattern1, Subpattern2, Subpattern3] = NSOJNMF_comodule(X1,X2,X3,W,H1,H2,H3,tt0,tt1,tt2,tt3);


module_file = 'Co_module_SOJNMF';
if ~isdir(module_file)
    mkdir(module_file);
end
cd('./Co_module_SOJNMF')

Co_mRNA = {};
for imRNA=1:K
    Co_mRNA(imRNA,1:length(Co_module{imRNA,2})) = mRNAs(1,Co_module{imRNA,2}+1);
end
xlswrite('Co_mRNA.xlsx', Co_mRNA);

Co_lncRNA = {};
for iLncRNA=1:K
    Co_lncRNA(iLncRNA,1:length(Co_module{iLncRNA,3})) = lncRNAs(1,Co_module{iLncRNA,3}+1);
end
xlswrite('Co_lncRNA.xlsx', Co_lncRNA);


Co_miRNA = {};
for imiRNA=1:K
    Co_miRNA(imiRNA,1:length(Co_module{imiRNA,4})) = miRNAs(1,Co_module{imiRNA,4}+1);
end
xlswrite('Co_miRNA.xlsx', Co_miRNA);

% save matrices
clearvars -except Co_module Subpattern1 Subpattern2 Subpattern3 OX1 OX2 OX3

save NSOJNMF_2234Comodule.mat


