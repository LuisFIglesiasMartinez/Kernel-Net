function [New_Models, L] = Add_a_variable_Numeric(Accepted_Models, Chosen_Models)

New_Models = zeros(1,size(Accepted_Models,2)+1);

for i = 1:size(Accepted_Models)

    Extra_Vars = ismember(Chosen_Models,Accepted_Models(i,:));
    
    Extra_Vars = find(Extra_Vars == 0);
    
    N = [repmat(Accepted_Models(i,:),[size(Extra_Vars,1),1]) Extra_Vars];
    
    A = [New_Models; N];
    
    New_Models = A;
    
end

New_Models = New_Models(2:end,:);

New_Models = unique((sort(New_Models')'), 'rows');

L = size(New_Models,1);
