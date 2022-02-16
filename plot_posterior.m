function [mode,mu,sigma,M,xdist,ydist] = plot_posterior(rounded_samples,xlab)
% FUNCTION PLOT_POSTERIOR
%
% plot inferred posterior distributions
% rounded_samples   === accepted samples rounded in day^-1 or day
% xlab              === x label


%%% PLOT FITTED HISTOGRAM
figure
hold on
h = histfit(rounded_samples,13,'kernel');
h(1).FaceColor = [245 245 245]./255;

%%% EXTRACT SUMMARY STATISTICS
[~,idx] = max(h(1).YData);
mode = h(1).XData(idx);
mu = mean(rounded_samples);
sigma = std(rounded_samples);
M = median(rounded_samples);

yL = get(gca,'YLim'); % yaxis for plotting
xdist = h(2).XData; % continous x variable
ydist = h(2).YData; % continuous y variable

%%% PLOT SUMMARY STATISTICS 
plot([mu mu],yL,'b','LineWidth', 2)
plot([mu+sigma mu+sigma],yL,'b--','LineWidth', 1.2)
plot([mu-sigma mu-sigma],yL,'b--','LineWidth', 1.2)
plot([M M],yL,'g','LineWidth', 1.2)
plot([mode mode],yL,'m','LineWidth', 1.2)
legend('discrete','continuous','\mu','\sigma^{+}','\sigma^{-}','median','mode')

if strcmp(xlab, '\beta_0')
    xlabel(['$' xlab '~(\mathrm{day}^{-1}$)'],'Interpreter', 'latex')
    title('People-people interactions')
elseif strcmp(xlab,'\gamma^{-1}')
    xlabel(['$' xlab '~(\mathrm{day}$)'],'Interpreter', 'latex')
    title('Incubation period')
elseif strcmp(xlab,'\delta^{-1}')
    xlabel(['$' xlab '~(\mathrm{day}$)'],'Interpreter', 'latex')
    title('Infectious period')
else
    xlabel(['$' xlab '~(\mathrm{virus}^{-1}$)'],'Interpreter', 'latex')
    title('People-pathogen interactions')
end

set(gca,'FontName', 'Times New Roman')  % Set it to times
set(gca,'FontSize', 16)
box on



end