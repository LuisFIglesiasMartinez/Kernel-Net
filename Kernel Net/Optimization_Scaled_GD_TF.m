function [SSE, Active_Fraction] = Optimization_Scaled_GD_TF(TF,mRNA_Targets)
%% Preprocesing

n = size(TF,1);

p = size(mRNA_Targets,2);

Y = mRNA_Targets.^-1;

B = zeros(3,p);

SSE = zeros(10000,1);

%% Estimate Starting Values
% Z is Active_Fraction^-1

Z = ones(n,1);

B2 = min(Y);

Y_1 = Y - repmat(B2,n,1);

B1 = sum(Y_1./repmat(TF,1,p))/(sum(TF.^2));

B(1,:) = B1';

B(2,:) = B2;

%% Do a Round of Steepest Descent

Search = 1;

e = 1e-4;

X = [Z./TF Z ones(n,1)];

    for j = 1:p 
         
    
         SSE(1) = SSE(1) + sum((Y(:,j) - X*B(:,j)).^2);
   
    end   
     
    t = 1;
    
    a = zeros(15,1);
    
    for i = 1:15
       
        a(i) = 1/(i^i);
    end

while Search
    
    t = t+1;
    %% Calculate the First Partial derivative of the Bs
    X = [Z./(TF) Z ones(n,1)];

    X_X = X'*X;

    dB = -2*X'*Y + 2*X_X*B;

    dB_D = dB;
    
    %% Scale it according to the Second Partial Derivative

    for q = 1:3
    
        dB_D(q,:) = dB(q,:)/(X_X(q,q));
  
    end

    %% Calculate the First Derivative of the Inverse of the Fraction 
    
    Y_1 = Y - repmat(B(3,:),n,1);

    W = (repmat(B(1,:),n,1) + repmat(B(2,:),n,1).*repmat(TF,1,p))./(repmat(TF,1,p));

    W_W = sum(W.*W,2);

    D = (2*W_W).^-1;

    dF = -2*sum(Y_1.*W,2)+2*W_W.*Z;

    %% Scale it 
    
    dF_D = dF.*D;
       
    SSE_1 = zeros(15,1);

    for i = 1:15
       
        Z_1 = Z - a(i)*dF_D;

        B_1 = B - a(i)*dB_D;
     
        B_1(1,B_1(1,:)<0) = 0;
    
        B_1(1,B_1(2,:)<0) = 0;
     
        Z_1(Z_1<1) = 1;
     
        X_1 = [Z_1./TF Z_1 ones(n,1)];
    
        for j = 1:p 
                  
            SSE_1(i) = SSE_1(i) + sum((Y(:,j) - X_1*B_1(:,j)).^2);
   
        end
                
    end
        [SSE(t),m] = min(SSE_1);    

    Z_1 = Z - a(m)*dF_D;

    B_1 = B - a(m)*dB_D;
    
     B_1(1,B_1(1,:)<0) = 0;
    
     B_1(1,B_1(2,:)<0) = 0;
     
     Z_1(Z_1<1) = 1;
     
     Z = Z_1;
     
     B = B_1;
     
     
        X = [Z_1./TF Z_1 ones(n,1)];

    if (SSE(t-1)-SSE(t))<e || t == 10000
                    
            B_LS = pinv(X'*X)*X'*Y;        
        
            
            
     B_LS(1,B_LS(1,:)<0) = 0;
    
     

     B_LS(1,B_LS(2,:)<0) = 0;
            
            
            SSE_1 = 0;
                             
            %% Do a Least Squares Estimation of the Bs
                  
            for j = 1:p 
         
         
                SSE_1 = SSE_1 + sum((Y(:,j) - X*B_LS(:,j)).^2);
            
            end   
            
            if SSE_1 < SSE(t)
                
                B = B_LS;
            
                SSE(t) = SSE_1;
            
            else
                
                 Search = 0;    
            end
            
    end
    
    if t == 10000 && Search
        
        SSE_Total = SSE;
        
        SSE = zeros(1000,1);
        
        SSE(1) = SSE_Total(end);
        
        t = 2;
        
    end
    
end


SSE = SSE(1:t);

Active_Fraction = Z.^-1;


