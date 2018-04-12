function Variable_Inclusion_Probability = Model_Search(Y,PCA,Posteriors,I,n,p,g,Prior,Var)

Variable_Inclusion_Probability = Posteriors;

Models = eye(p,p);

[Chosen_Models, ~] = Jeffreys_Model_Selection(Posteriors);

Log_Posteriors = Posteriors;

Best = max(Posteriors);

[P,Ranking] = sort(Posteriors, 'descend');

New_Best = Best;

Max_Model_Length = length(Chosen_Models);

if Max_Model_Length > 1
   
      Active = 1;  
            
      [New_Models, L] =  Make_Combinations(Chosen_Models,2);
      
    
while Active
            
        Posteriors = zeros(L,1);
        
        Accepted = zeros(L,1);
        
        for i = 1:L
        
            Vars = New_Models(i,:);
   
            X = Vector_Combination(PCA,Vars,Ranking,I);
            
            RSS = Var-Y'*X*pinv(X'*X)*X'*Y;
            
            k = size(X,2);
            
            Indexes = zeros(p,1);
            
            Indexes(Vars) = 1;
            
            not_i = find(Indexes == 0); 
            
            P = sum(log(Prior(Vars)))+sum(log(1-Prior(not_i)));
    
         Posteriors(i) = ((k+1)/2)*log(g/(g+1))+log(1/(g+1)*(RSS+g*Var))*(-((n-1)/2))+P;
    
            if Posteriors(i) > Best
                
            Accepted(i) = 1;
            
                if Posteriors(i)> New_Best
                    
               
                    New_Best = Posteriors(i);
            end
                
            end
            
        end
        
        if sum(Accepted) > 0
        
            Indexes = find(Accepted == 1);
            
        Posteriors = Posteriors(Indexes,:);
        
        New_Models = New_Models(Indexes,:);
            
         N_Models = zeros(size(New_Models,1),p);
        
        for i = 1:size(New_Models,1)
           
            
            N_Models(i,New_Models(i,:)) = 1;
            
        end
        
        M = [Models; N_Models];
        
        Models = M;        
        
        Lo  = [Log_Posteriors; Posteriors];
               
        Log_Posteriors = Lo;
        
        
           [Chosen_Models_1, ~] = Jeffreys_Model_Selection(Posteriors);
       
        
        New_Models = New_Models(Chosen_Models_1,:);
        
        
        
        [New_Models, L] = Add_a_variable_Numeric(New_Models, Chosen_Models); 
        
        
        Best = New_Best;
        
        else
            
            Active = 0;
            
        end
            
end
         
end

Posteriors = Logarithmic_Model_Averaging(Log_Posteriors);


Posteriors = repmat(Posteriors,1,p);
if size(Models,1) > 1
    
    
Variable_Inclusion_Probability = sum(Posteriors.*Models);

else
    
    Variable_Inclusion_Probability = Models;
    
end


