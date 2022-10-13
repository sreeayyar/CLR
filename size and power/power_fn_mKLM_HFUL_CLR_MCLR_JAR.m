function [p_mKLM, p_thli, p_CLR, p_MCLR, p_JAK]=power_fn_mKLM_HFUL_CLR_MCLR_JAR(n,m1,m2,k,rho,delta2,beta1,x)

grid_points = 2*x*5;

% Initialise p-value functions:
p_mKLM=zeros(grid_points,1);
p_thli=zeros(grid_points,1);
p_CLR=zeros(grid_points,1);
p_MCLR=zeros(grid_points,1);
p_JAK=zeros(grid_points,1);

% Generate the grid:
beta_lower_bound= beta1 - x;
beta_upper_bound= beta1 + x;
betagrid=beta_lower_bound(1):((beta_upper_bound(1)-beta_lower_bound(1))/(grid_points-1)):beta_upper_bound(1);


for i=1:grid_points
    [p_mKLM(i), p_thli(i), p_CLR(i), p_MCLR(i), p_JAK(i)]=sizes_mKLM_HFUL_CLR_MCLR_JAR(n,m1,m2,k,rho,delta2,betagrid(i),beta1);
end
