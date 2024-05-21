clc
clear all

% % Generate the Data

m       = 10;                                % Number of Agents
gamma   = 0.1;                                % Smoothing Parameter
mu      = 8*rand(1,m)-4;                      % Means of the local Gaussians  
sigma   = 0.1 + 0.5*rand(1,m);                % STD of the Local Gaussians

% load('new_30.mat')

x       = [-5:0.1:5];                         % Support of the discrete distribution
n       = size(x,2);                          % Cardinality of the support

% Generate the graph
typeG   = 'bi_path';    % Type of graph
G       = gen_graph(m,typeG,4/m) ;
bar_W   = full(laplacian(graph(G-diag(diag(G)))));
[~,D_w] = svd(bar_W); 
eigs_w = diag(D_w);
s_max = max(eigs_w);
bar_W = sparse(bar_W);
L       = s_max/gamma;                        % Lipschizt constant of the Gradients


N       = 5e3;
M       = 10;   %Samples from distributions
M2      = 10;   %Samples from communications

% load('F_Data_10.mat')
% load('F_Graph_10.mat')

L  = s_max/gamma;  
s_max = max(eigs_w);

%% Initializations
alpha(1) = 0;
beta(1)=0;
A(1) =0;

z = zeros(n,m);
lambda = zeros(n,m);
eta = zeros(n,m);
chi = zeros(n,m);
f = zeros(n,m);
aux = zeros(n,m);

cc = 3/2;
bb=1;
sigma_N=1e-1;
R=1;

bars = sqrtm(full(bar_W));

CCN  = zeros(n,m);
for i=1:N
    
    alpha(i+1) = (i+2)/2^cc;
    A(i+1) = A(i) + alpha(i+1);
    
    tau(i)=alpha(i+1)/A(i+1);
    
    beta(i+1)= L + (sigma_N/(2^(1/4)*sqrt(3)*R))*(i+2)^(cc);
    
    %% Compute sampled gradient and communication
    ran = randn(M,m);
    ran = bsxfun(@times, ran, sigma);
    ran = bsxfun(@plus, ran, mu);
    
    D = (1/6.5)*pdist2(reshape(ran,M*m,1),x','Euclidean').^2';
    CC = exp((kron(lambda,ones(1,M)) - D)/gamma);
    CC = CC./sum(CC);
    
    CC = reshape(CC',M,m*n);
    CC = mean(CC,1);
    CC = reshape(CC,m,n)';
    
    CCN  = sparse(zeros(n,m));
    for j=1:m
        indx = discretesample(CC(:,j), M2);
        for l=1:M2
            CCN(indx(l),j)=CCN(indx(l),j)+1;
        end
    end
    CCN = CCN/M2;
    
   
    %%
    
    aux = alpha(i)*CCN*bar_W +aux;
    z = -(1/beta(i+1))*aux;
%     z = -(1/beta(i+1))*(alpha(i+1)*CC*bar_W) - (beta(i)/beta(i+1))*z;
    lambda = tau(i)*z+(1-tau(i))*eta;
    
    %% Compute sampled gradient and communicatio
    ran = randn(M,m);
    ran = bsxfun(@times, ran, sigma);
    ran = bsxfun(@plus, ran, mu);
    
    D = (1/6.5)*pdist2(reshape(ran,M*m,1),x','Euclidean').^2';
    CC = exp((kron(lambda,ones(1,M)) - D)/gamma);
    CC = CC./sum(CC);
    
    CC = reshape(CC',M,m*n);
    CC = mean(CC,1);
    CC = reshape(CC,m,n)';
    
    CCN  = sparse(zeros(n,m));
    for j=1:m
        indx = discretesample(CC(:,j), M2);
        for l=1:M2
            CCN(indx(l),j)=CCN(indx(l),j)+1;
        end
    end
    CCN = CCN/M2;
    
    %% Continue
    chi = z -(alpha(i+1)/beta(i+1))*CCN*bar_W;
    eta = tau(i)*chi+(1-tau(i))*eta;
    
    f =  (alpha(i+1)*CCN + A(i)*f)/A(i+1);
    
     tot(i) = wass_smooth(gamma,x,lambda,mu,sigma,m,20);
     cons(i) = norm(f*bar_W)^2;
     
                 if mod(i,10)==0
%                 fifi{i} = f;
                plot(x,f)
%                         subplot(211)
%                      semilogx(tot)
%                      subplot(212)
%                     loglog(cons)
                drawnow;
            end
     
%      plot(tot)
%          drawnow;
    
%     if mod(i,10)==0
%         plot(f)
%         drawnow;
%     end

end