function beta_patch(ylims)
    % Add patches
    patch([13 20 20 13],[ylims(1) ylims(1) ylims(2) ylims(2)],'b','EdgeColor','None','FaceAlpha',.2);
    hold on
    patch([20 30 30 20],[ylims(1) ylims(1) ylims(2) ylims(2)],'r','EdgeColor','None','FaceAlpha',.2);
    hold on
end