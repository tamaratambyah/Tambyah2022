function plot_compare_data(x,xdata,ydata)
% FUNCTION PLOT_COMPARE_DATA
%
% x     === solution [S E I R C D]
% tspan === time array 
% xdata === data time array
% ydata === case data array


%%% PLOT S,E,I,R and I_fit
figure;
hold on
plot(xdata,x(:,1:4),'LineWidth', 1.2)
ax = gca; ax.ColorOrderIndex = 3;
plot(xdata,ydata,'o')
xlabel('$t$','Interpreter', 'latex')
ylabel('$N$','Interpreter', 'latex')
title('Fitted solution','Interpreter', 'latex')
legendInfo = {'$S$','$E$','$I$','$R$','$I^*$'};
legend(legendInfo,'Interpreter', 'latex')
set(gca,'FontName', 'Times New Roman')  % Set it to times
set(gca,'FontSize', 16)
box on

%%%% PLOT I and I_fit
figure;
hold on
ax = gca;
ax.ColorOrderIndex = 3;
plot(xdata,ydata,'o')
ax.ColorOrderIndex = 3;
plot(xdata,x(:,3),'LineWidth', 1.2)
xlabel('$t$','Interpreter', 'latex')
ylabel('$N$','Interpreter', 'latex')
title('Fitted solution','Interpreter', 'latex')
legendInfo = {'$I^*$','$I$'};
legend(legendInfo,'Interpreter', 'latex')
set(gca,'FontName', 'Times New Roman')  % Set it to times
set(gca,'FontSize', 16)
box on





end