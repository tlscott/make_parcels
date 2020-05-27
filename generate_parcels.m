function generate_parcels(FILENAMES,options)

%% generate parcels

% needed to run this:
% freesurfer installed and on your path, or set your path below
% an additional function called 'spm_ss_watershed.m' written by Alfonso
% Nieto-Castanon - now included within this script - 5/22/2020
% Code to create a sphere of voxels from Sung-Joo Lim.

% inputs - 
% FILENAMES = cell array where each cell contains a string indicating the full path to a subject's z-statistic file. 
% options = structure of the form below

% options = struct( ...
% 	'SET_FREESURFER', 'export FREESURFER_HOME=/Applications/freesurfer', ...
% 	'path_to_freesurfer_bin', ':/Applications/freesurfer/bin', ...
% 	'THRESH', 2.326348, ...
% 	'PCT_SUBJ_IN_PARC', 0.80, ...
%     'peak_spacing', 0, ...
% 	'EXPERIMENT', 'Nonword_Repetition', ...
% 	'PATH_TO_RESULTS_DIR', ['./' options.EXPERIMENT]);

% outputs:
% probability map (not scaled)
% probability map thresholded at 2 subjects 
% probability map thresholded at 2 subjects and then highly smoothed
% volume containing all local maxima in the smoothed map ('peaks')
% volume containing all parcels grown around local maxima
% volume containing only parcels with 80% or more subjects ('sig')
% additional useful data - how many subjs per parcel, how many voxels in
% each parcel, all input parameters, binary vector saying is a parcel is
% significant or not... etc.


%% Set things

% needed to find freesurfer functions particularly "mri_convert"; can skip if your path is already
% setup. you need freesurfer to run this code
% options.SET_FREESURFER = 'export FREESURFER_HOME=/Applications/freesurfer';
SET_FREESURFER = options.SET_FREESURFER;

% needed to access MATLAB freesurfer functions MRIread.m and MRIwrite.m. Replace with the location of those functions on your machine.
% options.path_to_freesurfer_bin = ':/Applications/freesurfer/bin';
path1 = getenv('PATH');
path1 = [path1 options.path_to_freesurfer_bin];
setenv('PATH', path1)

% z thresholds for reference: 
%   p = 0.01 : z = 2.326348
%   p = 0.001 : z = 3.090232
%   p = 0.0001 : z = 3.719015
%   p = 0.00001 : z = 4.264895

% options.Z_THRESH = 2.326348;
Z_THRESH = options.THRESH;

NSUBJ_THRESH = 2;

% options.PCT_SUBJ_IN_PARC = 0.80; % What percentage of subjects should contribute to a significant parcel?
PCT_SUBJ_IN_PARC = options.PCT_SUBJ_IN_PARC;

% options.EXPERIMENT = 'Nonword_Rep_NW4-NW1'; % just used to name output files like this "EXPERIMENT_probability_map.nii.gz", "EXPERIMENT_probability_map_thresh2subjs.nii.gz", etc.
EXPERIMENT = options.EXPERIMENT;

%% Set results directory 

% options.PATH_TO_RESULTS_DIR = './parcel_foldername';
RESULTS_DIR = [options.PATH_TO_RESULTS_DIR '_' date];
mkdir(RESULTS_DIR)

%% Open up new blank volume

new_volume = MRIread([FILENAMES{1}],0);
zeros_vol = zeros(size(new_volume.vol));
new_volume.vol = zeros_vol;

%% Binarize maps

for i = 1:length(FILENAMES)
    
    dummy_vol = zeros_vol;
    ind_subj_vol = MRIread([FILENAMES{i}],0);
    sig_voxels = find(ind_subj_vol.vol >= Z_THRESH);
    dummy_vol(sig_voxels) = 1;
    
    new_volume.vol = new_volume.vol + dummy_vol;
    
    clear ind_subj_vol sig_voxels dummy_vol
    
end

MRIwrite(new_volume,[RESULTS_DIR '/' EXPERIMENT '_probability_map.nii.gz'],'float')

%% Threshold at 2 subjects

new_volume.vol(new_volume.vol < NSUBJ_THRESH) = 0;

MRIwrite(new_volume,[RESULTS_DIR '/' EXPERIMENT '_probability_map_thresh' num2str(NSUBJ_THRESH) 'subjs.nii.gz'],'float')

%% Smooth thresholded volume

FREESURFER_COMMAND = ['mri_convert -i ' RESULTS_DIR '/' EXPERIMENT '_probability_map_thresh' num2str(NSUBJ_THRESH) 'subjs.nii.gz -o ' ...
    RESULTS_DIR '/' EXPERIMENT '_probability_map_thresh' num2str(NSUBJ_THRESH) 'subjs_smoothed.nii.gz --fwhm 6'];
unix([SET_FREESURFER ' && ' FREESURFER_COMMAND]);

% Temporary:
smoothed_vol = MRIread([RESULTS_DIR '/' EXPERIMENT '_probability_map_thresh' num2str(NSUBJ_THRESH) 'subjs_smoothed.nii.gz'],0);

%% Watershed algorithm

IDX = find(smoothed_vol.vol > NSUBJ_THRESH);

[D,P] = spm_ss_watershed(-smoothed_vol.vol,IDX,options.peak_spacing);

parcel_vol = new_volume;
peak_vol = new_volume;

parcel_vol.vol = D;
peak_vol.vol = P;

MRIwrite(parcel_vol,[RESULTS_DIR '/' EXPERIMENT '_probability_map_thresh' num2str(NSUBJ_THRESH) 'subjs_smoothed_parcels.nii.gz'],'float');
MRIwrite(peak_vol,[RESULTS_DIR '/' EXPERIMENT '_probability_map_thresh' num2str(NSUBJ_THRESH) 'subjs_smoothed_peaks.nii.gz'],'float');

%% How many subjects have significant voxels in each parcel?

[~,~,parcel_numbers] = find(unique(D));

subjs_in_parcel_mat = zeros(length(FILENAMES),length(parcel_numbers));

for i = 1:length(FILENAMES)
    
    dummy_vol = zeros_vol;
    ind_subj_vol = MRIread([FILENAMES{i}],0);
    sig_voxels = find(ind_subj_vol.vol >= Z_THRESH);
    dummy_vol(sig_voxels) = 1;
    
    for j = 1:length(parcel_numbers)
       
        temp = find((dummy_vol > 0) & (parcel_vol.vol == parcel_numbers(j)));
        
        if length(temp) > 0
            subjs_in_parcel_mat(i,j) = 1;
        end
        
        clear temp
        
    end
    
    clear ind_subj_vol sig_voxels dummy_vol
    
end

count_in_parcel = sum(subjs_in_parcel_mat,1);

for i = 1:length(parcel_numbers)    
    nvoxels_in_parcel(i) = length(find(parcel_vol.vol == parcel_numbers(i)));   
end

sig_parcels = find(count_in_parcel >= PCT_SUBJ_IN_PARC*length(FILENAMES));

is_significant = zeros(size(nvoxels_in_parcel));
is_significant(sig_parcels) = 1;

sig_parcel_volume = D;
sig_parcel_volume(~ismember(D,sig_parcels)) = 0;

sig_parcel = parcel_vol;
sig_parcel.vol = sig_parcel_volume;

MRIwrite(sig_parcel,[RESULTS_DIR '/' EXPERIMENT '_probability_map_thresh' num2str(NSUBJ_THRESH) 'subjs_smoothed_parcels_sig.nii.gz'],'float');

[~,peak_labels,max_intersect_labels] = labelParcels(P,D,parcel_numbers);

is_significant = is_significant';
count_in_parcel = count_in_parcel';
nvoxels_in_parcel = nvoxels_in_parcel';

parcel_report = table(parcel_numbers,is_significant,count_in_parcel,nvoxels_in_parcel,peak_labels,max_intersect_labels);

writetable(parcel_report,[RESULTS_DIR '/' EXPERIMENT '_parcel_report.csv'])

save([RESULTS_DIR '/' EXPERIMENT '_addtnl_data.mat']);

end

function [D,P] = spm_ss_watershed(A,IDX,rad_n_vox)
% SPM_SS_WATERSHED watershed segmentation
%
% C=spm_ss_watershed(A);
% C=spm_ss_watershed(A,idx);
%

% note: assumes continuous volume data (this implementation does not work well with discrete data). In practice this means having sufficiently-smoothed volume data
%

% New! Check for other peaks within a sphere of a certain radius
[offsets,~] = construct_sphere(rad_n_vox);

sA=size(A);

%zero-pad&sort
if nargin<2, IDX=find(~isnan(A)); IDX=IDX(:); else IDX=IDX(:); end
[~,idx]=sort(A(IDX)); idx=IDX(idx); 
[pidx{1:numel(sA)}]=ind2sub(sA,idx(:));
pidx=mat2cell(1+cat(2,pidx{:}),numel(pidx{1}),ones(1,numel(sA)));
eidx=sub2ind(sA+2,pidx{:});
sA=sA+2;
N=numel(eidx);

%neighbours (max-connected; i.e. 26-connected for 3d)
[dd{1:numel(sA)}]=ndgrid(1:3);
d=sub2ind(sA,dd{:});
d=d-d((numel(d)+1)/2);d(~d)=[];

%assigns labels
C=zeros(sA);P=zeros(sA-2);
m=1;
for n1=1:N
    c=C(eidx(n1)+d);
    c=c(c>0);
    temp_sphere_idxs = calc_sphere_idxs(A,idx(n1),offsets);
    peak_sphere_idxs=find(P(temp_sphere_idxs),1);
    if isempty(c) && isempty(peak_sphere_idxs)
        C(eidx(n1))=m;P(idx(n1))=m;m=m+1;
    elseif isempty(c) && ~isempty(peak_sphere_idxs)
        C(eidx(n1))=P(temp_sphere_idxs(peak_sphere_idxs));
    elseif ~any(diff(c))
        C(eidx(n1))=c(1);
    end
    
end
D=zeros(size(A));D(idx)=C(eidx);

end



function [offsets,distances] = construct_sphere(radius_fix)

% construct a sphere in a voxel space (isotrophic size).
% input: radius of the sphere - # of voxels
% 
% output: sphere coordinates in voxelspace with origin of 0,0,0 

    radius = radius_fix;


    side=ceil(radius)*2+1;

    % make an array that is sufficienly large (actually too big)
    offsets=zeros(side^3,3);

    % keep track of where we store indices in offsets:
    % for every location within the sphere (in the code below),
    % add one to row_pos and then store the locations
    row_pos=0;

    % the grid positions (relative to the origin) to consider
    single_dimension_candidates=floor(-radius):ceil(radius);

    % Consider the candidates in all three spatial dimensions using
    % nested for loops and an if statement.
    % For each position:
    % - see if it is at most at distance 'radius' from the origin
    % - if that is the case, increase row_pos and store the position in
    %   'offsets'
    %
    % >@@>
    for x=single_dimension_candidates
        for y=single_dimension_candidates
            for z=single_dimension_candidates
                if x^2+y^2+z^2<=radius^2
                    row_pos=row_pos+1;
                    offsets(row_pos,:)=[x y z];
                end
            end
        end
    end
    % <@@<

    % cut off empty values at the end
    offsets=offsets(1:row_pos,:);

    % compute distances
    unsorted_distances=sqrt(sum(offsets.^2,2));

    % sort distances and apply to offsets
    [distances,idxs]=sort(unsorted_distances);
    offsets=offsets(idxs,:);
end

function sphere_idxs = calc_sphere_idxs(vol1,temp_center,offsets)

[x,y,z] = ind2sub(size(vol1),temp_center);
    
    %% Step 3: Draw a sphere around it
    
    sphere_idxs = [];
    n = 1;
    
    for i = 1:size(offsets,1)
        
        if ((x+offsets(i,1)) <= size(vol1,1)) && ((y+offsets(i,2)) <= size(vol1,2)) && ((z+offsets(i,3)) <= size(vol1,3)) ...
                && ((x+offsets(i,1)) > 0) && ((y+offsets(i,2)) > 0) && ((z+offsets(i,3)) > 0)
            sphere_idxs(n) = sub2ind(size(vol1),x+offsets(i,1),y+offsets(i,2),z+offsets(i,3));
            n = n+1;
        end
        
    end
end