function Network = Kernel_Net(mRNA_Exp,TFs, Prior,Sigma)

%% Function to Infer Gene Regulatory Networks from Gene Expression Measurements

%% Input Description

% mRNA_Exp:  Gene expression matrix (PxN) N is the number of genes and P is
% the number of patients.

% Optional 

% TFs: Vector (1xK)Indexes of the genes which are transcription factors
% K.It should match the order in which the genes in mRNA_Exp matrix are.

% Prior: an NxP matrix with the prior probabilities that each transcription
% factor regulates a target gene. The columns are transcription factors and
% the rows are target genes. 

% Sigma: a hyperparameter for the RBF Kernel

%% Output Description

% Network: an NxP matrix with the posterior probabilities that each
% interaction takes place. 


%% Preprocessing Steps 

N = size(mRNA_Exp); 

P = N(1); % No. of Patients.

N = N(2); % No. of Genes in the set.

mRNA_Exp = mRNA_Exp - repmat(mean(mRNA_Exp),N,1);

mRNA_Exp = mRNA_Exp./repmat(std(mRNA_Exp),N,1);


switch nargin
    
    case 1
        
        K = size(mRNA_Exp,2);
        
        TFs = 1:K;
        
        Prior = ones(N,K)*0.5;
        
        Sigma = 100;
        
    case 2
        
  
        K = length(TFs); % No. of Transcription Factors.
       
        Prior = ones(N,K)*0.5;
        
        Sigma = 100;
        
    case 3
        
        Sigma = 100;
        
end
        
          
     

% 1 Scaling the Gene Expression Data and calculating some useful variables



Network = zeros(N,K); % Preallocate memory for the Output Network

% 2 Calculating the Kernel matrix between transcription factor gene
% expression and the Kernel Principal Components for each of the transcription
% factors

g = 1/N;

PCA = cell(K,1);

Posteriors = zeros(N,K);

V = sum(mRNA_Exp.*mRNA_Exp);

for i = 1:K
    
  
%% Get the Squared Distances

X = mRNA_Exp(:,TFs(i));

D2 = distSquared(X,X);

Kernel = exp(-D2/Sigma);
    
  
   PCA{i,1} = Kernel_PCA(Kernel);
   
  % Making the Principal Components Orthogonal
   
   PCA{i,1} =  PCA{i,1} - repmat(mean(PCA{i,1}),P, 1);

   RSS = sum((mRNA_Exp-PCA{i,1}*pinv(PCA{i,1}'*PCA{i,1})*PCA{i,1}'*mRNA_Exp).^2);

    k = size(PCA{i,1},2);
    
    if i == 1
        
        not_i = 2:K;
        
    elseif i < K
        
        not_i = [1:i-1 i+1:K];
        
    else
        
        not_i = 1:K-1;
        
    end
    
    Prior_K = log(Prior(:,i)')+sum(log(1-Prior(:,not_i)'));
    
    Posteriors(:,i) = (k/2)*log(g/(g+1))+log((RSS'+g*V')/(1+g))*(-((P-1)/2))+Prior_K';

    Posteriors(TFs(i),i) = -Inf;

   
end

%% Variable Selection Steps
Is = cell(N,1);

Net = cell(N,1);

for i = 1:N
    
    Is{i} = find(~isinf(Posteriors(i,:)));
    
    Net{i} = Model_Search(mRNA_Exp(:,i),PCA,Posteriors(i,Is{i})',Is{i},P,length(Is{i}),g,Prior(i,Is{i}),V(i));
    
end

for i = 1:N 
    
    
   Network(i, Is{i}) = Net{i};
    
end

