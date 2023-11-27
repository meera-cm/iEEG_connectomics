function [niix, mapping, dummy] = get_niix(filenames,Dt)

    for i = 1:size(filenames,1)
        try
            nii = ea_load_nii(filenames{i});

            if i == 1
                % allocate space to speed up the process
                niix = zeros(size(filenames,1),size(nii.img(:),1));
                mapping = zeros(Dt.nchannels,1);
            end

            niix(i,:) = nii.img(:);
            [~,fname]=fileparts(nii.fname);

            % julian probs also wrote a function for this..
            cname = strfind(fname,'cname'); 
            pos = strfind(fname,'_pos'); 
            name = fname(cname(1)+6:pos(1)-1);

            % add the betavalues for that filename (so for that channel)
            ichname = find(ismember(Dt.info.ChannelName, name));
            mapping(i) = ichname;
        catch
            disp('error');
        end
    end

    % somehow theres Inf values in niix
    niix(niix == Inf) = NaN;
    dummy = nii;
    dummy.fname = 'dummy';
end