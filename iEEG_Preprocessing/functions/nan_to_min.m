function vals = nan_to_min(vals)
    for i=1:size(vals,2)
        for j=1:size(vals,3)
            vals(isnan(vals(:,i,j)),i,j) = min(vals(:,i,j));
        end
    end
end