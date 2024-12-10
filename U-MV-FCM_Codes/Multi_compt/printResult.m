function [AR] = printResult(X, label, K, kmeansFlag)

if kmeansFlag == 1
    indic = litekmeans(X, K, 'Replicates',20);
else
    [~, indic] = max(X, [] ,2);
end
result = bestMap(label, indic);
[ac, nmi_value, cnt] = CalcMetrics(label, indic);

[AR,RI]=RandIndex(label, indic);
Outs =valid_external(label, indic);
AR=1-ErrorRate(label, indic,K)/size(indic,1)

NMI  =nmi(label, indic)
JI=Outs(3) 
FMI=Outs(4)

disp(sprintf('AR: %0.4f\t%d/%d\tRI:%0.4f\t', AR, cnt, length(label), RI));
