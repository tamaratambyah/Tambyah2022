function plot_accepted_samples(rounded_samples,parms_fit,errors)
% FUNCTION PLOT_ACCEPTED_SAMPLES
%
% plot accepted \beta_0,\gamma^{-1},\delta^{-1} samples as scatter plot
% samples_accepted  === accepted samples from ABC rounded to days
% parms_fit         === parameter set with best fit 
% errors            === array of discrepancy values
% num_accepted      === number of accepted samples from ABC

figure
hold on
scatter3(rounded_samples(:,1),1./rounded_samples(:,2),1./rounded_samples(:,3),...
    50,errors(:),'filled')
scatter3(parms_fit(1),parms_fit(2),parms_fit(3),60,'r','filled')
view(3)
grid on
colorbar
title(['Accepted samples: ' num2str(size(rounded_samples,1))])
xlabel('$\beta_0$','Interpreter', 'latex')
ylabel('$\gamma^{-1}$','Interpreter', 'latex')
zlabel('$\delta^{-1}$','Interpreter', 'latex')
set(gca,'FontName', 'Times New Roman')  
set(gca,'FontSize', 16)
box on

end