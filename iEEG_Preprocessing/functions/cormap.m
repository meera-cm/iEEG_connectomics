function cormap(niix,betavalue,dummy,range_name,outfolder,xtype)
    
    if ~exist('xtype','var')
        xtype = '';
    end
    
    dummy.fname     = sprintf('%s/%s/corr_%s.nii', outfolder, ...
        range_name, xtype);

    if isempty(dir(dummy.fname))
        [~,~,rall,~,~,~] = wjn_fox_loom(niix,betavalue,'pearson',0);
        dummy.img(:)    = rall;
        ea_write_nii(dummy);
    end
    
end