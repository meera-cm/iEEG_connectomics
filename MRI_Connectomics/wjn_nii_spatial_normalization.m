
function nii = wjn_nii_spatial_normalization(nii)

if ischar(nii) || iscell(nii)
    nii = ea_load_nii(nii);
end

nii = ea_load_nii(nii);
nii.fname = ['normalized_' nii.fname];
nii.img(:)=nii.img(:)./nanmean(nii.img(:));
ea_write_nii(nii)



%  serotonin_nii = ea_load_nii('02_5HT1b_az_hc36_beliveau.nii')
%  serotonin_nii.fname = ['normalized_' serotonin_nii.fname];
%  serotonin_nii.img(:)=serotonin_nii.img(:)./nanmean(serotonin_nii.img(:));
%  ea_write_nii(serotonin_nii)