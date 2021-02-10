function relativeAbundances = simulateretronfitness_single(parameters,relativeAbundances)
 
% This function expects two structures as inputs
% parameters:
%   parameters.r is the editing rate per generation
%   parameters.f is the fitness of the edited strain
%   parameters.editingLength duration in generations that editing is
%       maintained
%   parameters.totalLength duration in generation of the full simulation
 
% relativeAbundances:
%   relativeAbundances.WT is the fraction of WT strains
%   relativeAbundances.unedited is the fraction of unedited strains
%   relativeAbundances.edited is the fraction of edited strains
%
%could initialize with an actual number of cells... this would tell you
%when a given mut would become extinct
for iGeneration=1:parameters.totalLength

    relativeAbundances.Bene_edit(end+1,1)=relativeAbundances.Bene_edit(end) ...
        +relativeAbundances.Bene_unedit(end) * parameters.r;
    relativeAbundances.Bene_unedit(end+1,1)=relativeAbundances.Bene_unedit(end) ...
        -relativeAbundances.Bene_unedit(end) * parameters.r;


    relativeAbundances.Neut(end+1,1)=relativeAbundances.Neut(end);

        
    %update frequencies based on fitness
    
    %note: this assumes no clonal interference, that is, relative abundance is
    %purely a function of growth rate, and cells never approach carrying
    %capacity nor compete
    
    %This is the change in abundance of the edited population relative to the
    %WT/unedited population which is assumed to have a fitness of 1
    relativeAbundances.Bene_edit(end)=relativeAbundances.Bene_edit(end) * parameters.BeneF; %grow with rate f
    relativeAbundances.Bene_unedit(end)=relativeAbundances.Bene_unedit(end) * 1; %neutral fitness
   
    relativeAbundances.Neut(end)=relativeAbundances.Neut(end) * 1; %Neut has neutral fitness
    

    
end
    



