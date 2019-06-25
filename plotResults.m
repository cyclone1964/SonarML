% This script plots the results of the investigation into various
% SVM models
[Kernels,Params,Costs,Errors] = textread('results.dat','%s%f%f%f','delimiter','\t');

useNamedFigure('linear'); clf
Indices = find(strcmp(Kernels,'linear'));
plot(Costs(Indices),Errors(Indices),'o-');
set(gca,'YLim',[0 2 ]);
set(gca,'XScale','log');
xlabel('Cost (log)');
ylabel('CVE (%)');
title('Linear Models');
prettyPlot;
print('-dpng','LinearErrors.png');
useNamedFigure('polynomial'); clf
for Degree = [3 5 7]
    Indices = find(strcmp(Kernels,'polynomial') & Params == Degree);
    hold on;
    plot(Costs(Indices),Errors(Indices),'o-')
end
set(gca,'XScale','log');
set(gca,'YLim',[0 20]);
xlabel('Cost');
ylabel('CVE (%)')
legend('Order 3','Order 5','Order 7');
title('Polynomial CVE');
prettyPlot
print('-dpng','PolynomialErrors.png');
useNamedFigure('radial'); clf;

for Gamma = [0.001 0.01 0.1 1]
    Indices = find(strcmp(Kernels,'radial') & Params == Gamma);
    hold on;
    plot(Costs(Indices),Errors(Indices),'o-');
end
set(gca,'XScale','log');xlabel('Cost');
set(gca,'YLim',[0 5]);
ylabel('CVE (%)')
legend('Gamma 0.001','Gamma 0.01','Gamma 0.1','Gamma 1');
title('Radial CVE');
prettyPlot
print('-dpng','RadialErrors.png')




