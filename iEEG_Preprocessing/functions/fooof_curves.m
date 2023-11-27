function [L, G] = fooof_curves(Dt)
    %% compute periodic and aperiodic part from fooof parameters
    L = zeros(Dt.nchannels,Dt.nfrequencies);
    G = zeros(Dt.nchannels,Dt.nfrequencies);

    for i=1:Dt.nchannels

        c = Dt.info.peaks.gaussian_params(Dt.info.peaks.ic==i,1);
        a = Dt.info.peaks.gaussian_params(Dt.info.peaks.ic==i,2);
        w = Dt.info.peaks.gaussian_params(Dt.info.peaks.ic==i,3);
        b = Dt.info.peaks.background_params(i,1);

        if size(Dt.info.peaks.background_params,2) == 2
            k = 0;
            x = Dt.info.peaks.background_params(i,2);
        else
            k = Dt.info.peaks.background_params(i,2);
            x = Dt.info.peaks.background_params(i,3);
        end

        L(i,:) = b - log10(k + Dt.frequencies.^x);

        for j = 1:size(c,1)
            G(i,:) = G(i,:) + a(j)*gaussmf(Dt.frequencies, [w(j) c(j)])';
        end
    end
end