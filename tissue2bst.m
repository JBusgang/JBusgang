% Tissue Probability Map Import to 5-Tissue for BST
% Written by Jacob Busgang and Edouard Delaire

function varargin = tissue2bst(segIn, subjectname, protocol_anatdir, Origin, extract2compare)

% INPUT: 
%   - segIn            = segmentation file of format .nii.gz or brainstorm 
%                        anatomy to be converted to 5 tissue segmentation.
%   - subjectname      = Subject name in bst database
%   - protocol_anatdir = Full directory to the protocol's anat folder.
%   - Origin           = describes segIn's map type, with options:
%                        'CHARM' : CHARM 10-tissue segmentation.
%                        'ASEG'  : ASEG 44-tissue segmentation.
%                        'BST'   : Brainstorm 5-tissue segmentation.
%    - extract2compare  = optional input 0 or 1. 1 tells the function to 
%                         create varargin as a workspace variable copy of 
%                         the .Cube data in order to compare to other maps.
% OUTPUT:
%    - varargin         = extracted Cube data in 5-tissue 256x256x256
%                         single array format.
%    - segmentation_5tissues_(Origin) is created upon using the function
%                                     to import a non-BST file.
%
% This function is designed to work with the Brainstorm software
% https://neuroimage.usc.edu/brainstorm

% Set segmentation labels to the 5 expected bst labels
Labels = {[0],'Background',[0 0 0];[1],'White',[220 220 220];[2],'Gray',[130 130 130];[3],'CSF',[44 152 254];[4],'Skull',[255 255 255];[5],'Scalp',[255 205 184]};

% Check extract2compare
x2c = exist('extract2compare');

%% CHARM to Brainstorm 5 tissues transform
if strcmp(Origin, 'CHARM')
    sVol = in_mri(segIn, 'Nifti1'); % Load segmentation
    cubeNew = cast(sVol.Cube, 'single');
    for i = 1:256
        for j = 1:256
            for k = 1:256
                if cubeNew(i,j,k)==6 || cubeNew(i,j,k)==9 ||cubeNew(i,j,k)==10
                    cubeNew(i,j,k)=5; % Eyeballs, blood, muscle to Scalp. Circle back on making blood CSF, or fifth type.
                elseif cubeNew(i,j,k)==7 || cubeNew(i,j,k)==8
                    cubeNew(i,j,k)=4; % Compact bone and spongy bone to Skull
                end
            end
        end    
    end

    sVol.Cube = cast(cubeNew, 'uint8');
    sVol.Labels = Labels;
    if x2c && extract2compare
        varargin = cubeNew;
    end
end

%% ASEG to Brainstorm 5 tissues transform
if strcmp(Origin, 'ASEG')
    sVol = in_mri(segIn, 'BST'); % Load segmentation
    cubeNew = cast(sVol.Cube, 'single')
    for i = 1:256
        for j = 1:256
            for k = 1:256
                % ismember if statements are cleaner, but very slow compared to multiple or's.
                if cubeNew(i,j,k)==3 || cubeNew(i,j,k)==8 || cubeNew(i,j,k)==10 || cubeNew(i,j,k)==11 || cubeNew(i,j,k)==12 ...
                    || cubeNew(i,j,k)==16 || cubeNew(i,j,k)==17 || cubeNew(i,j,k)==18 || cubeNew(i,j,k)==26 || cubeNew(i,j,k)==42 ... 
                    || cubeNew(i,j,k)==47 || cubeNew(i,j,k)==49 || cubeNew(i,j,k)==50 || cubeNew(i,j,k)==51 || cubeNew(i,j,k)==52 ...
                    || cubeNew(i,j,k)==53 || cubeNew(i,j,k)==54 || cubeNew(i,j,k)==58
                    cubeNew(i,j,k)=2; 
                elseif cubeNew(i,j,k)==2 || cubeNew(i,j,k)==7 || cubeNew(i,j,k)==28 || cubeNew(i,j,k)==31 || cubeNew(i,j,k)==41 ...
                        || cubeNew(i,j,k)==46 || cubeNew(i,j,k)==60 || cubeNew(i,j,k)==63 || cubeNew(i,j,k)==77 || cubeNew(i,j,k)==85 ...
                        || cubeNew(i,j,k)==251 || cubeNew(i,j,k)==252 || cubeNew(i,j,k)==253 || cubeNew(i,j,k)==254 || cubeNew(i,j,k)==255
                    cubeNew(i,j,k)=1;
                elseif cubeNew(i,j,k)==4 || cubeNew(i,j,k)==5 || cubeNew(i,j,k)==13 || cubeNew(i,j,k)==14 || cubeNew(i,j,k)==15 ...
                        || cubeNew(i,j,k)==24 || cubeNew(i,j,k)==30 || cubeNew(i,j,k)==43 || cubeNew(i,j,k)==44 || cubeNew(i,j,k)==62
                    cubeNew(i,j,k)=3;
                else
                    cubeNew(i,j,k)=0;
                end 
            end
        end
    end

    sVol.Cube = cast(cubeNew, 'uint8');
    sVol.Labels = Labels;
    if x2c && extract2compare
        varargin = cubeNew;
    end
end

% Import Brainstorm 5 tissues segmentation file
if strcmp(Origin,'BST')
    sVol = in_mri(segIn, 'BST'); % Load segmentation
    cubeNew = cast(sVol.Cube, 'single');
    if x2c && extract2compare
        varargin = cubeNew;
    end
end

% Save to brainstorm
if ~strcmp(Origin, 'BST')
    [sSubject, iSubject] = bst_get('Subject', subjectname);
    c = clock;
    strTime = sprintf('_%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
    vol_fn = fullfile(protocol_anatdir, subjectname, ['subjectimage_T1_segmentation_', Origin, strTime,'_volatlas.mat']);


    vol_name = file_unique(sprintf('segmenation_5tissues_%s', Origin), {sSubject.Anatomy.Comment}, 1);

    add_vol_data(sVol,...
                 cubeNew , ...
                 vol_fn , ...
                 vol_name, ...
                iSubject, ...
                sprintf('Segmentation from %s to 5 tissue', Origin));
end
end
%



function sSubject = add_vol_data(sMri, data, vol_fn, vol_comment, iSubject, history_comment)
%% Save a volume map to brainstorm with given data

ProtocolSubjects = bst_get('ProtocolSubjects');

if (iSubject == 0) % Default subject
	sSubject = ProtocolSubjects.DefaultSubject;
else % Normal subject
    sSubject = ProtocolSubjects.Subject(iSubject);
    if sSubject.UseDefaultAnat
        sSubject = ProtocolSubjects.DefaultSubject;
        iSubject = 0;
    end 
end

if ~exist(vol_fn,'file')
    is_new = 1;
else
    is_new = 0;
end

sMri.Cube       = data;
sMri.Comment    = vol_comment;
sMri.Histogram  = mri_histogram(sMri.Cube);


if nargin > 5
    sMri = bst_history('add', sMri, 'import', history_comment);
end

% Save new MRI in Brainstorm format
out_mri_bst(sMri, vol_fn,'v6');

% If the MRI already exist, no need to go further.
if ~is_new 
    return;
end
%% ===== STORE NEW MRI IN DATABASE ======
% New anatomy structure
iAnatomy = find(cell2mat(cellfun(@(a) ~isempty(a), ...
                   strfind({sSubject.Anatomy.FileName}, vol_fn),...
                   'uni', false )));
% Add new anatomy
if isempty(iAnatomy)
    iAnatomy = length(sSubject.Anatomy) + 1;
end
sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
sSubject.Anatomy(iAnatomy).FileName = file_short(vol_fn);
sSubject.Anatomy(iAnatomy).Comment  = vol_comment;

% == Update database ==
% Default subject
if (iSubject == 0)
	ProtocolSubjects.DefaultSubject = sSubject;
% Normal subject 
else
    ProtocolSubjects.Subject(iSubject) = sSubject;
end
bst_set('ProtocolSubjects', ProtocolSubjects);

%% ===== UPDATE GUI =====
% Refresh tree
panel_protocols('UpdateNode', 'Subject', iSubject);
panel_protocols('SelectNode', [], 'subject', iSubject, -1 );
% Save database
db_save();
end
