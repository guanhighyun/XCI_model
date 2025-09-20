% The time point at which the average amount of bound Xist over 10
% simulations are calculated
sample_time_point = [40, 100, 190, 1000];

for i = 1:numel(sample_time_point)
    plot_average_bound_Xist(sample_time_point(i))
end

function plot_average_bound_Xist(sample_time_point)
load('Fig4D.mat','Xist_count','SPEN_count');
figure;

for i = 1:10
    curr_Xist_count = Xist_count{i}; % current realization
    Xist_count_at_final_time_point = curr_Xist_count(end,:); 

    % Calculate total Xist count over all time points from X1 (binding site
    % 1-100), X2 (binding site 101-200), and X3 (binding site 201-300).
    sum_Xist_count = [sum(Xist_count_at_final_time_point(1:100)),...
        sum(Xist_count_at_final_time_point(101:200)),...
        sum(Xist_count_at_final_time_point(201:300))];

    % Find the chromosome with maximum total Xist count. 
    % This is an inactive X chromosome covered by a lot of bound Xist and SPEN.
    [~,max_idx] = max(sum_Xist_count);

    % Extract the bound Xist ditribution from the inactive X.
    cell_Xist_count = {curr_Xist_count(sample_time_point,1:100),curr_Xist_count(sample_time_point,101:200),curr_Xist_count(sample_time_point,201:300)};
    sum_max_Xist_count(i,:) = (cell_Xist_count{max_idx});
end

plot([60,60],[1.1,1.2],'r','LineWidth',5);hold on;
bar(1:100,mean(sum_max_Xist_count),0.1,'LineWidth',2); 
hold off; ylim([0,1.2]); xlim([0,100]); set(gca,'fontsize',25)
end
