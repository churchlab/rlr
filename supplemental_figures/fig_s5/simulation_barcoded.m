close all
clear all

%Parameter valuations
legendHolder={};
figure,
for i= [.001 .01 .05 .07 .09 .11] %unknown
    parameters.r=i;
    %edited to only need one value of f
    parameters.f=0;
    parameters.editingLength=40;%experiment variable
    parameters.totalLength=40;%experiment variable
    %initial conditions
    relativeAbundancesbar.WT=1.0;%experiment variable
    relativeAbundancesbar.edited=0;%experiment variable
    relativeAbundancesbar.unedited=1.0;%experiment variable
    relativeAbundancesbar.barcoded=relativeAbundancesbar.edited+relativeAbundancesbar.unedited; %this is to visualize total barcoded set
    relativeAbundancesbar.barcoded_over_total = relativeAbundancesbar.barcoded(end) / ...
            (relativeAbundancesbar.barcoded(end) + relativeAbundancesbar.WT(end));
    
    relativeAbundancesbar=simulateretronfitness_barcodes(parameters,relativeAbundancesbar);
    
    plot(relativeAbundancesbar.barcoded_over_total);
    
    
    legendHolder{end+1}=num2str(i);
    hold all
end

legend(legendHolder);
title(['Barcode abundance where f is ',num2str(parameters.f)]);
yL = get(gca,'YLim');
%line([parameters.editingLength parameters.editingLength],yL,'Color','r');
xlabel('Generations'); % x-axis label
ylabel('Relative Abundance'); % y-axis label

        
