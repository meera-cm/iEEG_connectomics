function plot_stat_cor(vals,keys,ranges)
    %% plot peak attribute distributions
    figure('units','normalized','outerposition',[0 0 1 1]);
    s=1;
    for n=1:size(vals,3)
        for i=1:size(vals,3)
            subplot(size(vals,3),size(vals,3),s);
            s=s+1;
            cors = [];
            
            for j=1:size(vals,2)
                scatter(vals(:,j,i),vals(:,j,n));
                cors = [cors, corrcoef(vals(:,j,i), ...
                    vals(:,j,n), 'rows','complete')];
                hold on;
            end
         
            if n == size(vals,3)
                xlabel(keys(i));
            end
            
            if i == 1
                ylabel(keys(n));
            end
            title(strcat('R = ', sprintf(' %.2f,',cors(1,2:2:end))));
            %legend(ranges);
        end
    end
    
end