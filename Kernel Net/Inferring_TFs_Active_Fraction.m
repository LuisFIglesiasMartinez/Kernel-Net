function [TF_Activity, Active_Set] = Inferring_TFs_Active_Fraction(Net,Data,TF_Index,T)

Data = Data - min(min(Data)) + 0.01;


switch nargin
    case 2
        K = size(Data,2);
        TF_Index = 1:K;

        T = 0.1232;

    case 3
        T = 0.1232;
        
end

TF = Data(:,TF_Index);


TF_Activity = TF;

Active_Set = zeros(1,length(TF_Index));

for i = 1: length(TF_Index)
    
    Ones = find(Net(:,i)>T);
    
    if length(Ones)>0
    %% Calculate Correlation
    
    Active_Set(i) = 1;
    
     mRNA_Targets = Data(:,Ones);
    
     p = zeros(length(Ones),1);
    
     for j = 1:length(Ones)
  
         p(j) = corr(TF(:,i),mRNA_Targets(:,j));
    
     end
    
     Pos = find(p>0);
     
     Negs = find(p<0);
     
    if (length(Pos)*length(Negs))>0 
        
        disp('first if')
    
    [~,TF_Activity(:,i)] = Optimization_Scaled_GD_TF_Inh_Act(TF(:,i),mRNA_Targets);
    
    elseif length(Pos)>0
        
        
        [~,TF_Activity(:,i)] = Optimization_Scaled_GD_TF(TF(:,i),mRNA_Targets);
        
    else
         [~,TF_Activity(:,i)] = Optimization_Scaled_GD_TF_Inh_2(TF(:,i),mRNA_Targets);
        
    
    end
    
    end
    
    
    
end

TF_Activity = TF_Activity(:,Active_Set);

Active_Set = find(Active_Set ==1);

