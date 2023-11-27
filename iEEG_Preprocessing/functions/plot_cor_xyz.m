function plot_cor_xyz(vals,keys,ranges,coords)
    %% plot peak attribute distributions
    figure('units','normalized','outerposition',[0 0 1 1]);
    s=1;
    cnames = [string('x') string('y') string('z')];
    for n=1:size(coords,2)
        for i=1:size(vals,3)
            subplot(size(coords,2),size(vals,3),s);
            s=s+1;
            cors = [];
            
            for j=1:size(vals,2)
                scatter(coords(:,n),vals(:,j,i));
                cors = [cors, corrcoef(coords(:,n), ...
                    vals(:,j,i), 'rows','complete')];
                hold on;
            end
         
            if n == size(coords,2)
                xlabel(keys(i));
            end
            
            if i == 1
                ylabel(cnames(n));
            end
            title(strcat('R = ', sprintf(' %.2f,',cors(1,2:2:end))));
            %legend(ranges);
        end
    end
    
end