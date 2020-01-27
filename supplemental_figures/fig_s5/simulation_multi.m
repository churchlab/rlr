

%%%% lol will have to move this outside the loop for this case
parameters.BeneF = 1.02;
parameters.DetriF = 0.9;
parameters.LethalF = 0;
parameters.r = 0.05;

parameters.editingLength=20; %experiment variable
parameters.totalLength=40; %experiment variable
parameters.totalpopsize=1;
%initial conditions

relativeAbundances.Neut=1;%experiment variable
relativeAbundances.Bene_edit=0.00;%experiment variable
relativeAbundances.Bene_unedit=1;%experiment variable
relativeAbundances.Detri_edit=0.00;%experiment variable
relativeAbundances.Detri_unedit=1;%experiment variable
relativeAbundances.Lethal_edit=0.00;%experiment variable
relativeAbundances.Lethal_unedit=1;%experiment variable

%specify initial barcode abundance
barcodeAbundances.Neut = relativeAbundances.Neut;
barcodeAbundances.Bene = relativeAbundances.Bene_edit + relativeAbundances.Bene_unedit;
barcodeAbundances.Detri = relativeAbundances.Detri_edit + relativeAbundances.Detri_unedit;
barcodeAbundances.Lethal = relativeAbundances.Lethal_edit + relativeAbundances.Lethal_unedit;

%run simulation
relativeAbundances=simulateretronfitness_multi(parameters,relativeAbundances);

%this depicts the total barcoded set for a given mutation
barcodeAbundances.Neut=relativeAbundances.Neut;
barcodeAbundances.Bene=relativeAbundances.Bene_edit+ relativeAbundances.Bene_unedit;
barcodeAbundances.Detri=relativeAbundances.Detri_edit+ relativeAbundances.Detri_unedit;
barcodeAbundances.Lethal=relativeAbundances.Lethal_edit+ relativeAbundances.Lethal_unedit;

barcodeFrequencies.Neut = barcodeAbundances.Neut/barcodeAbundances.Neut;
barcodeFrequencies.Bene = barcodeAbundances.Bene/barcodeAbundances.Neut;
barcodeFrequencies.Detri = barcodeAbundances.Detri/barcodeAbundances.Neut;
barcodeFrequencies.Lethal = barcodeAbundances.Lethal/barcodeAbundances.Neut;


%plots
legendHolder={};
figure,
for iStrain=fieldnames(barcodeFrequencies)'
    iStrain=cell2mat(iStrain);
    plot(barcodeFrequencies.(iStrain));
    legendHolder{end+1}=iStrain;
    hold all
end

legend(legendHolder);
title(['r is ',num2str(parameters.r),' and f is',num2str(parameters.f)]);
yL = get(gca,'YLim');
line([parameters.editingLength
 parameters.editingLength],yL,'Color','r')
        

xlabel('Generations'); % x-axis label
ylabel('Relative Abundance'); % y-axis label