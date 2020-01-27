function relativeAbundancesbar = simulateretronfitness(parameters,relativeAbundancesbar)
 
% This function expects two structures as inputs
% parameters:
%   parameters.r is the editing rate per generation
%   parameters.f is the fitness of the edited strain
%   parameters.editingLength duration in generations that editing is
%       maintained
%   parameters.totalLength duration in generation of the full simulation
 
% relativeAbundancesbar:
%   relativeAbundancesbar.WT is the fraction of WT strains
%   relativeAbundancesbar.unedited is the fraction of unedited strains
%   relativeAbundancesbar.edited is the fraction of edited strains
%
%could initialize with an actual number of cells... this would tell you
%when a given mut would become extinct
for iGeneration=1:parameters.totalLength
    if parameters.editingLength>0
        %Update frequences based on editing %%%%%%%%this happens when gen <
        %editingLength
        parameters.editingLength=parameters.editingLength-1; %a generation happens
        
        relativeAbundancesbar.edited(end+1,1)=relativeAbundancesbar.edited(end) ...
            +relativeAbundancesbar.unedited(end) * parameters.r;
        
        relativeAbundancesbar.unedited(end+1,1)=relativeAbundancesbar.unedited(end) ...
            -relativeAbundancesbar.unedited(end) * parameters.r;
    else
        relativeAbundancesbar.edited(end+1,1)=relativeAbundancesbar.edited(end);
        relativeAbundancesbar.unedited(end+1,1)=relativeAbundancesbar.unedited(end);
    end
        relativeAbundancesbar.WT(end+1,1)=relativeAbundancesbar.WT(end);
        
    %update frequencies based on fitness
    
    %note: this assumes no clonal interference, that is, relative abundance is
    %purely a function of growth rate, and cells never approach carrying
    %capacity nor compete
    
    %This is the change in abundance of the edited population relative to the
    %WT/unedited population which is assumed to have a fitness of 1
    relativeAbundancesbar.edited(end)=relativeAbundancesbar.edited(end)+relativeAbundancesbar.edited(end) * parameters.f;
    relativeAbundancesbar.unedited(end)=2*relativeAbundancesbar.unedited(end);
    relativeAbundancesbar.WT(end)=2*relativeAbundancesbar.WT(end);
    
    %This is the change in the total abundance of the entire population
    totalPopulationAbundance = relativeAbundancesbar.edited(end) ...
        + relativeAbundancesbar.unedited(end) + relativeAbundancesbar.WT(end);
    
    %This nornalizes the population back to its orginal size
    relativeAbundancesbar.edited(end)=relativeAbundancesbar.edited(end)/(totalPopulationAbundance/2);
    relativeAbundancesbar.unedited(end)=relativeAbundancesbar.unedited(end)/(totalPopulationAbundance/2);
    relativeAbundancesbar.WT(end)=relativeAbundancesbar.WT(end)/(totalPopulationAbundance/2);
    
    %this depicts the total barcoded set for a given mutation
    relativeAbundancesbar.barcoded(end+1,1)=relativeAbundancesbar.edited(end)+relativeAbundancesbar.unedited(end);
    
    %this depicts the fraction barcoded, relative to total reads
    relativeAbundancesbar.barcoded_over_total(end+1,1) =  relativeAbundancesbar.barcoded(end) / ...
            (relativeAbundancesbar.barcoded(end) + relativeAbundancesbar.WT(end));
end
    



