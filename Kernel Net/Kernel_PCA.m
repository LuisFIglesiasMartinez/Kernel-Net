function Projections = Kernel_PCA(Kernel)

%% Step 1 Eigendecomposition

[Vectors, Values] = eig(Kernel);

%% Step 2 Normalisation 

cutoff = 0.01;

Values = diag(Values);

Indx = find(Values > cutoff);

Indx = sort(Indx, 'descend');

Values = Values(Indx);

Vectors = Vectors(:,Indx);

Vectors = Vectors ./repmat(sum(Vectors.^2),size(Vectors,1),1);

Vectors = Vectors.*((repmat(Values', size(Vectors,1),1).^-(1/2)));

% Step 3: Projections

Projections = Kernel*Vectors;
