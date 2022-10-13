%% Power Graphs
% This code uses the power_fn_mKLM_HFUL_CLR_MCLR_JAR function to generate
% the four power graphs in the paper. Note, only power of MCLR, JAR and
% mKLM are included in the graphs, although the code also generates power 
% curves for CLR and H-LIML. 

% set seed
rng(1);

%how big do you want grid to be:
x=5;

%True beta:
beta1=0;

%% Graph 1: power_n200_k30_rho02_d10_mKLM_CLR_MCLR_JAR

[p_mKLM, p_thli, p_CLR, p_MCLR, p_JAK]=power_fn_mKLM_HFUL_CLR_MCLR_JAR(200,5000,5000,30,0.2,10,beta1,x);

% Save figure:
beta1=0;
% Generate the grid:
grid_points= 2*x*5; 
beta_lower_bound= beta1 - x;
beta_upper_bound= beta1 + x;
betagrid=beta_lower_bound(1):((beta_upper_bound(1)-beta_lower_bound(1))/(grid_points-1)):beta_upper_bound(1);

figurename = strcat('power_n200_k30_rho02_d10_mKLM_CLR_MCLR_JAK', '.eps');
fig = figure(1);
betagrid_plot = betagrid;
plot(-betagrid_plot',p_mKLM(ismember(betagrid,betagrid_plot)),'-.',...
    -betagrid_plot,p_MCLR(ismember(betagrid,betagrid_plot))', '-.',...
    -betagrid_plot,p_JAK(ismember(betagrid,betagrid_plot))', '-.',...
    'LineWidth',2);

legend({sprintf('mKLM'),sprintf('MCLR'),sprintf('Jackknife AR')},...
    'Location','southoutside','Orientation','horizontal','FontSize',14 );
legend('boxoff');
xlabel('\beta'); xlim([min(betagrid_plot) max(betagrid_plot)])
ylabel('Rejection Probability of 5% Test, H_0:\beta=0','FontSize',14 );
ylim([0 1])
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,figurename,'-depsc','-r0')

%% Graph 2: power_n200_k30_rho02_d2_mKLM_CLR_MCLR_JAR

[p_mKLM, p_thli, p_CLR, p_MCLR, p_JAK]=power_fn_mKLM_HFUL_CLR_MCLR_JAR(200,5000,5000,60,0.2,30,beta1,x);

% Save figure:
beta1=0;
% Generate the grid:
grid_points= 2*x*5;
beta_lower_bound= beta1 - x;
beta_upper_bound= beta1 + x;
betagrid=beta_lower_bound(1):((beta_upper_bound(1)-beta_lower_bound(1))/(grid_points-1)):beta_upper_bound(1);

figurename = strcat('output/power_n200_k60_rho02_d30_mKLM_CLR_MCLR_JAK', '.eps');
fig = figure(1);
betagrid_plot = betagrid;
plot(-betagrid_plot',p_mKLM(ismember(betagrid,betagrid_plot)),'-.',...
    -betagrid_plot,p_MCLR(ismember(betagrid,betagrid_plot))', '-.',...
    -betagrid_plot,p_JAK(ismember(betagrid,betagrid_plot))', '-.',...
    'LineWidth',2);

legend({sprintf('mKLM'),sprintf('MCLR'),sprintf('Jackknife AR')},...
    'Location','southoutside','Orientation','horizontal','FontSize',14 );
legend('boxoff');
xlabel('\beta'); xlim([min(betagrid_plot) max(betagrid_plot)])
ylabel('Rejection Probability of 5% Test, H_0:\beta=0','FontSize',14 );
ylim([0 1])
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,figurename,'-depsc','-r0')

%% Graph 3: power_n200_k30_rho02_d30_mKLM_CLR_MCLR_JAR

[p_mKLM, p_thli, p_CLR, p_MCLR, p_JAK]=power_fn_mKLM_HFUL_CLR_MCLR_JAR(200,5000,2500,30,0.2,30,beta1,x);

% Save figure:
beta1=0;
% Generate the grid:
grid_points= 2*x*5;
beta_lower_bound= beta1 - x;
beta_upper_bound= beta1 + x;
betagrid=beta_lower_bound(1):((beta_upper_bound(1)-beta_lower_bound(1))/(grid_points-1)):beta_upper_bound(1);
 
figurename = strcat('output/power_n200_k30_rho02_d30_mKLM_CLR_MCLR_JAK', '.eps');
fig = figure(1);
betagrid_plot = betagrid;
plot(-betagrid_plot',p_mKLM(ismember(betagrid,betagrid_plot)),'-.',...
    -betagrid_plot,p_MCLR(ismember(betagrid,betagrid_plot))', '-.',...
    -betagrid_plot,p_JAK(ismember(betagrid,betagrid_plot))', '-.',...
    'LineWidth',2);

legend({sprintf('mKLM'),sprintf('MCLR'),sprintf('Jackknife AR')},...
    'Location','southoutside','Orientation','horizontal','FontSize',14 );
legend('boxoff');
xlabel('\beta'); xlim([min(betagrid_plot) max(betagrid_plot)])
ylabel('Rejection Probability of 5% Test, H_0:\beta=0','FontSize',14 );
ylim([0 1])
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,figurename,'-depsc','-r0')

%% Graph 4: power_n200_k30_rho02_d2_mKLM_CLR_MCLR_JAR

[p_mKLM, p_thli, p_CLR, p_MCLR, p_JAK]=power_fn_mKLM_HFUL_CLR_MCLR_JAR(200,5000,5000,30,0.2,60,beta1,x);

% Save figure:
beta1=0;
% Generate the grid:
grid_points= 2*x*5;
beta_lower_bound= beta1 - x;
beta_upper_bound= beta1 + x;
betagrid=beta_lower_bound(1):((beta_upper_bound(1)-beta_lower_bound(1))/(grid_points-1)):beta_upper_bound(1);

figurename = strcat('output/power_n200_k30_rho02_d60_mKLM_CLR_MCLR_JAK', '.eps');
fig = figure(1);
betagrid_plot = betagrid;
plot(-betagrid_plot',p_mKLM(ismember(betagrid,betagrid_plot)),'-.',...
    -betagrid_plot,p_MCLR(ismember(betagrid,betagrid_plot))', '-.',...
    -betagrid_plot,p_JAK(ismember(betagrid,betagrid_plot))', '-.',...
    'LineWidth',2);

legend({sprintf('mKLM'),sprintf('MCLR'),sprintf('Jackknife AR')},...
    'Location','southoutside','Orientation','horizontal','FontSize',14 );
legend('boxoff');
xlabel('\beta'); xlim([min(betagrid_plot) max(betagrid_plot)])
ylabel('Rejection Probability of 5% Test, H_0:\beta=0','FontSize',14 );
ylim([0 1])
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,figurename,'-depsc','-r0')

