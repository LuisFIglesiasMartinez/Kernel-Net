function [Chosen_Models, Log_Posteriors] = Jeffreys_Model_Selection(M_Ls)

[Max] = max(M_Ls);

Chosen_Models = find(Max<= M_Ls+2.3);

Log_Posteriors = M_Ls(find(Max<=M_Ls+2.3));



