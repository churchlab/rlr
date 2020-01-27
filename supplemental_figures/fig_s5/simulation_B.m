%similation for "B" with multiple traces on one plot
legendHolder={};
figure,

%Parameter valuations
plot_set = {'edited'};


for i= [.05] %[.05 .07 .09 .11]
    parameters.r=i;
    for j=1 %: 0.2 : 1.0 %unknown
        parameters.f=j;
        parameters.editingLength=100;%experiment variable
        parameters.totalLength=100;%experiment variable
        parameters.totalpopsize=1;
        %initial conditions
        relativeAbundances.WT=0.0;%experiment variable
        relativeAbundances.edited=0.00;%experiment variable
        relativeAbundances.unedited=1.0;%experiment variable
        relativeAbundances.barcoded=relativeAbundances.edited +relativeAbundances.unedited ; %this is to visualize total barcoded set
 
    relativeAbundances=simulateretronfitness(parameters,relativeAbundances);
        
        %legendHolder={};
        %figure,
        %loop through strains and plot their curves
        
        %This loop plots all values of iStrain
        %for iStrain=fieldnames(relativeAbundances)'
        for iStrain=plot_set
            iStrain=cell2mat(iStrain);
            plot(relativeAbundances.(iStrain));
            %legendHolder{end+1}=iStrain; %this plots edited etc
            legendHolder{end+1}= num2str(i)
            hold all
        end

    end
end

legend(legendHolder);
title(['r is ',num2str(parameters.r),' and f is',num2str(parameters.f)]);
yL = get(gca,'YLim');
        

xlabel('Generations'); % x-axis label
ylabel('Relative Abundance'); % y-axis label