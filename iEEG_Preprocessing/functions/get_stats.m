function [vals, val_names] = get_stats(Dt,pp,pi,range)

    val_names = [string('number of peaks'), string('maximum peak amplitude'), ...
        string('integral of highest peak'), string('total integral')];
    
    %% get peak data for low and high beta
    for j=1:size(range,1)
        for i=1:Dt.nchannels

            pn = find(Dt.info.peaks.ic == i);

            il = find(Dt.info.peaks.peak_params(pn,1) > range(j,1) & ...
                Dt.info.peaks.peak_params(pn,1) <= range(j,2)) + min(pn) - 1;

            [maxl, argmaxl] = max(pp(il,2));

            if ~isempty(il)
                vals(i,j,:) = [size(il,1) maxl pi(il(argmaxl)) sum(pi(il))];
            else
                vals(i,j,:) = [0 NaN NaN NaN];
            end
        end
    end
end