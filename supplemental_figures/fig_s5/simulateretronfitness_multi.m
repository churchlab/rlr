function relativeAbundances = simulateretronfitness_multi(parameters,relativeAbundances)
 
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
    if parameters.editingLength>0 
        parameters.editingLength=parameters.editingLength-1; %a generation happens
        
        relativeAbundances.Bene_edit(end+1,1)=relativeAbundances.Bene_edit(end) ...
            +relativeAbundances.Bene_unedit(end) * parameters.r;
        relativeAbundances.Bene_unedit(end+1,1)=relativeAbundances.Bene_unedit(end) ...
            -relativeAbundances.Bene_unedit(end) * parameters.r;
        
        relativeAbundances.Detri_edit(end+1,1)=relativeAbundances.Detri_edit(end) ...
            +relativeAbundances.Detri_unedit(end) * parameters.r;
        relativeAbundances.Detri_unedit(end+1,1)=relativeAbundances.Detri_unedit(end) ...
            -relativeAbundances.Detri_unedit(end) * parameters.r;
        
        relativeAbundances.Lethal_edit(end+1,1)=relativeAbundances.Lethal_edit(end) ...
            +relativeAbundances.Lethal_unedit(end) * parameters.r;
        relativeAbundances.Lethal_unedit(end+1,1)=relativeAbundances.Lethal_unedit(end) ...
            -relativeAbundances.Lethal_unedit(end) * parameters.r;
        
    else
        relativeAbundances.Bene_edit(end+1,1)=relativeAbundances.Bene_edit(end);
        relativeAbundances.Bene_unedit(end+1,1)=relativeAbundances.Bene_unedit(end);
        
        relativeAbundances.Detri_edit(end+1,1)=relativeAbundances.Detri_edit(end);
        relativeAbundances.Detri_unedit(end+1,1)=relativeAbundances.Detri_unedit(end);
        
        relativeAbundances.Lethal_edit(end+1,1)=relativeAbundances.Lethal_edit(end);
        relativeAbundances.Lethal_unedit(end+1,1)=relativeAbundances.Lethal_unedit(end);
  
    end
        relativeAbundances.Neut(end+1,1)=relativeAbundances.Neut(end);

        
    %update frequencies based on fitness
    
    %note: this assumes no clonal interference, that is, relative abundance is
    %purely a function of growth rate, and cells never approach carrying
    %capacity nor compete
    
    %This is the change in abundance of the edited population relative to the
    %WT/unedited population which is assumed to have a fitness of 1
    relativeAbundances.Bene_edit(end)=relativeAbundances.Bene_edit(end) * parameters.BeneF; %grow with rate f
    relativeAbundances.Bene_unedit(end)=relativeAbundances.Bene_unedit(end) * 1; %neutral fitness
   
    relativeAbundances.Detri_edit(end)=relativeAbundances.Detri_edit(end) * parameters.DetriF; %grow with rate f
    relativeAbundances.Detri_unedit(end)=relativeAbundances.Detri_unedit(end) * 1; %neutral fitness
    
    relativeAbundances.Lethal_edit(end)=relativeAbundances.Lethal_edit(end) * parameters.LethalF; %grow with rate f
    relativeAbundances.Lethal_unedit(end)=relativeAbundances.Lethal_unedit(end) * 1; %neutral fitness
    
    relativeAbundances.Neut(end)=relativeAbundances.Neut(end) * 1; %Neut has neutral fitness
    
    %This is the change in the total abundance of the entire population
    %is this section still required if i'm multiplying by one? maybe helps
    %to eliminate drift due to rounding?
    %totalPopulationAbundance = relativeAbundances.Bene_edit(end) + relativeAbundances.Bene_unedit(end) ...
    %    + relativeAbundances.Detri_edit(end) + relativeAbundances.Detri_unedit(end) ...
    %    + relativeAbundances.Lethal_edit(end) + relativeAbundances.Lethal_unedit(end) ...
    %    + relativeAbundances.Neut(end);
    
    %This nornalizes the population back to its orginal size
    %relativeAbundances.edited(end)=relativeAbundances.edited(end)/(totalPopulationAbundance/parameters.totalpopsize);
    %relativeAbundances.unedited(end)=relativeAbundances.unedited(end)/(totalPopulationAbundance/parameters.totalpopsize);
    %relativeAbundances.WT(end)=relativeAbundances.WT(end)/(totalPopulationAbundance/parameters.totalpopsize);
    

    
    %these are the actual frequencies, normalized to Neut
    %totalabundance = 
    
end
    



