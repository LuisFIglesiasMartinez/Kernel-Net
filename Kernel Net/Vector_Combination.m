function Vectors = Vector_Combination(Vects,Vars,Ranking,I)


%% First find the Variable with the highest Ranking


Order =I(find(ismember(Ranking,Vars)==1));

Vectors = Vects{Order(1),1};

for i = 2:length(Order)

    V = Vects{Order(i),1};
    
    V = [Vectors V];
    
    Vectors = V;
   
end



