function  BMA = Logarithmic_Model_Averaging(Log_Posteriors)



Log_Rep = repmat(Log_Posteriors,1,size(Log_Posteriors,1));

Diff = Log_Rep-Log_Rep';

Fractions = exp(Diff);

Sum_of_Fractions = sum(Fractions);

BMA = (Sum_of_Fractions.^-1)';




