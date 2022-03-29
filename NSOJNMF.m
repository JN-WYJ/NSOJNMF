function [W,H1,H2,H3] = NSOJNMF(X1,X2,X3,theta1,theta2,theta3,R12,R13,R23,K,lambda1,lambda2,alpha,gamma,maxiter,speak,fid,iloop)
%
% Multiple NMF using euclidean distance update equations:
%
% INPUT:
% X1 (N,M1): N (dimensionallity) x M1 (samples) non negative input matrix
% X2 (N,M2): N (dimensionallity) x M2 (samples) non negative input matrix
% X3 (N,M3): N (dimensionallity) x M3 (samples) non negative input matrix
% theta1 (M1,M1): the interaction matrix of X1-X1
% theta2 (M2,M2): the interaction matrix of X2-X2
% theta3 (M3,M3): the interaction matrix of X3-X3
% R12 (M1,M2): the interaction matrix of X1-X2
% R13 (M1,M3): the interaction matrix of X1-X3
% R23 (M2,M3): the interaction matrix of X2-X3
% K        : Number of components
% alpha    : orthogonal parameter
% lambda1  : the weighted parameter of the network regularization constraint
% lambda2  : the weighted parameter of the network regularization constraint
% gamma   : sparse parameter
% maxiter  : Maximum number of iterations to run
% speak    : prints iteration count and changes in connectivity matrix 
%            elements unless speak is 0
%
% OUTPUT:
% W        : N x K matrix
% H1       : K x M1 matrix
% H2       : K x M2 matrix
% H3       : K x M3 matrix
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_iter = 20; % iterations between print on screen and convergence test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X1 X2 and X3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (min(min(X1)) < 0) || (min(min(X2)) < 0) || (min(min(X3)) < 0)
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for same rows in X1 X2 and X3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,m1] = size(X1);
[n2,m2] = size(X2);
[n3,m3] = size(X3);

if (n1~=n2) || (n1~=n3)   
    error('Input matrices should have the same rows');
    return
end
n = n1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W, H1, H2 and H3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = rand(n,K);
H1 = rand(K,m1);
H2 = rand(K,m2);
H3 = rand(K,m3);

I = eye(K);
E1 = ones(K,m1);
E2 = ones(K,m2);
E3 = ones(K,m3);


for iter = 1:maxiter
    % multiplicative update method
    H1 = H1.*(W'*X1 + 2*alpha*H1 + lambda1*H1*theta1 + 1/2*lambda2*(H2*R12' + H3*R13'))./((W'*W)*H1 + 2*alpha*H1*H1'*H1 + gamma/2*E1 + eps);     
    H2 = H2.*(W'*X2 + 2*alpha*H2 + lambda1*H2*theta2 + 1/2*lambda2*(H1*R12 + H3*R23'))./((W'*W)*H2 + 2*alpha*H2*H2'*H2 + gamma/2*E2 + eps);   
    H3 = H3.*(W'*X3 + 2*alpha*H3 + lambda1*H3*theta3 + 1/2*lambda2*(H1*R13 + H2*R23))./((W'*W)*H3 + 2*alpha*H3*H3'*H3 + gamma/2*E3 + eps);    
    W = W.*([H1 H2 H3]*[X1 X2 X3]')'./(W*([H1 H2 H3]*[H1 H2 H3]') + eps);   
  
    if (rem(iter,print_iter) == 0) & speak 
		% compute residues
		diff1 = nmf_euclidean_dist(X1,W*H1) + nmf_euclidean_dist(X2,W*H2) + nmf_euclidean_dist(X3,W*H3);
		diff2 = alpha*nmf_euclidean_dist(H1*H1',I) + alpha*nmf_euclidean_dist(H2*H2',I) + alpha*nmf_euclidean_dist(H3*H3',I);
		diff3 = (-lambda1*trace(H1*theta1*H1')) + (-lambda1*trace(H2*theta2*H2')) + (-lambda1*trace(H3*theta3*H3'));
		diff4 = (-lambda2*trace(H1*R12*H2')) + (-lambda2*trace(H1*R13*H3')) + (-lambda2*trace(H2*R23*H3'));
		diff5 = gamma*norm(H1,1) + gamma*norm(H2,1) + gamma*norm(H3,1);
		diff = diff1 + diff2 + diff3 + diff4 + diff5;
		errorx1 = mean(mean(abs(X1-W*H1)))/mean(mean(X1));
		errorx2 = mean(mean(abs(X2-W*H2)))/mean(mean(X2));
		errorx3 = mean(mean(abs(X3-W*H3)))/mean(mean(X3));
		errorx = errorx1 + errorx2 + errorx3;
    
		disp(['Iter = ',int2str(iter),...
				', relative error = ',num2str(errorx)])
	
		fprintf(fid,'%s\n',[sprintf('nloop = \t'),int2str(iloop),...
				sprintf('\t Iter = \t'),int2str(iter),...
				sprintf('\t diff = \t'),num2str(diff),...
				sprintf('\t relative error = \t'),num2str(errorx)]);
        
		if errorx < 10^(-5), break, end
	end
    
end

end

function err = nmf_euclidean_dist(X,Y)

err = sum(sum((X-Y).^2));
end
