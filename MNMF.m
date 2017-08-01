function [U, M, H, C, L] = MNMF(S, M, U, H, C, B1, B2, alpha, beta, lambda)

% The demo is written by Xiao Wang (wangxiao007@mail.tsinghua.edu.cn), and the details of the algorithm can be
% found in "Community Preserving Network Embedding" (AAAI 2017).

%%--------------input-----------------
% S: S1+5*S2, the first- and second-order proximities (n-by-n);
% M: the initialized basis matrix (n-by-m);
% U: the initialized representations of nodes (n-by-m);
% H: the initialized community indicator matrix (n-by-k);
% C: the initialized representations of communities (k-by-m);
% B1: the adjacency matrix B1(i,i)=0 (n-by-n);
% B2: (k_i*k_j)/2e (n-by-n);
% alpha, beta, lambda: the values of parameters, usually alpha and beta need to be tuned and we can set lambada = 1e9;
%%-------------output------------------
% U: the optimal representations of nodes;
% M: the optimal basis matrix;
% H: the optimal community indicator matrix;
% C: the optimal representations of communities;
% L: the final values of objective function.
%%-------------------------------------

K = size(H,2);
I = eye(K);
X = U';
for i = 1:200
    % update M:
    M = M.*((S*U)./max(realmin, M*(U'*U)));
    
    % update U:
    X = X.*((M'*S+alpha*C'*H')./max(realmin, (M'*M+alpha*(C'*C))*X));
    U = X';
    
    % update C:
    C = C.*((H'*U)./max(realmin, C*U'*U));
    
    % update H:
    B2H = B2*H;
    HHH = H*(H'*H);
    B1H = B1*H;
    UC = U*C';
    sqroot1 = sqrt((2*beta*B2H).^2+16*lambda*HHH.*(2*beta*B1H+2*alpha*UC+(4*lambda-2*alpha)*H));
    sqroot2 = sqrt((-2*beta*B2H+sqroot1)./max(realmin,(8*lambda*HHH)));
    H = H.*sqroot2;
end

% objective function values
first = norm(S-M*U','fro')^2;
second = alpha*norm(H-U*C','fro')^2;
third = -beta*trace(H'*(B1-B2)*H);
constraint = lambda*norm(H'*H-I,'fro')^2;
L = first + second + third + constraint;