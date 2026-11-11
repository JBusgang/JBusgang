function vargout = process_sensitivity_probmap(vargin)

%
% _________________________________________________________________________
% process_sensitivity_probmap
% Compute the probability of coverage across the cortex surface.
% _________________________________________________________________________
% INPUTS
%  - N summed sensitivites of group, coregistered to template, in form of
%    1xN cell. 
% _________________________________________________________________________
% Written by Jacob Busgang, Concordia University (2026)

sFiles = vargin

% check for matching size and 2 or more files
if length(sFiles)<2;
    error("Requires more than one sensitivity map.")
end
surf1 = in_bst_channel(sFiles{1}).SurfaceFile;
for i = 1:(length(sFiles))
    surf2 = in_bst_channel(sFiles{i}).SurfaceFile;
    if (surf1~=surf2)
        error("Surface files not the same.")
    end
end

% Binarize and sum sensitivities across maps
sTotalAmp = zeros(length(in_bst_channel(sFiles{1}).ImageGridAmp), 1);
for i = 1:length(sFiles)
    sAmp = in_bst_channel(sFiles{i}).ImageGridAmp;
    sAmp(sAmp~=0)=1;
    sTotalAmp = sTotalAmp + sAmp;
end

% Normalize for probabilistic model
sTotalAmp = sTotalAmp/max(sTotalAmp)

% Save as map
sOutput = in_bst_channel(sFiles{1});
sOutput.ImageGridAmp = sTotalAmp;
sOutput.Comment = 'Group Spatial Coverage Probability Map';
vargout = sOutput;