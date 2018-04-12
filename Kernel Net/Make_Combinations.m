function [New_Models, L] =  Make_Combinations(Variables,Level)

New_Models = combnk(Variables,Level);

L = size(New_Models,1);





