% REPLICATION: Mikusheva and Sun (2022) EMPIRICAL SIMULATION
% ADDITIONS: homoskedastic, linear data-generating process; add mKLM, CLR and MCLR
% for comparison. Focus on size of tests only. 
% Maintain structure of the replication code!

% This code produces a power comparison of M-CLR, J-AR and mKLM tests using
% a simulated dataset that is inspired by Angrist & Krueger (1991), and
% draws on the simulation design in Angrist & Fraseden (2022). 

clear

%critical values:
alpha = 0.05;
mKLM95=chi2inv(1-alpha,1);
JAK95=norminv(1-alpha,0,1);
t_hli95=norminv(1-alpha/2,0,1);

beta = 0.1; % ground truth

% Generate the grid:
beta_lower_bound= -5;
beta_upper_bound= 5;
betagrid=beta_lower_bound(1):0.25:beta_upper_bound(1);
grid_points = length(betagrid);

%load data:
load simulation_AK91_no_ZW_homolinear_tenth.dat; % run qobsims_reverse_engineer.do to generate this dataset
simulation_AK91 = simulation_AK91_no_ZW_homolinear_tenth;

%extract random sub-sample and keep it fixed
rng(1); 
random = rand(size(simulation_AK91,1),1);
sh = 0.0025; % randomly select a subsample that is 0.25% of the full sample size
select = random <= sh;

%% Create Instruments

simulation_AK91 = simulation_AK91(select == 1,:);
simulation_AK91 = sortrows(simulation_AK91,[5 6 7]);

control = simulation_AK91(:,3);
%sigz_normalized = simulation_AK91(:,2);
educbar = simulation_AK91(:,1);

qob = simulation_AK91(:,5);
yob = simulation_AK91(:,6);
pob = simulation_AK91(:,7);

obs = size(control,1);

%create state of birth and year of birth dummies:
SOB_index=[2 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 ...
    32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56]; %drop 56 to avoid multicollinartity
SOB_dummies=zeros(obs,length(SOB_index));
for n=1:length(SOB_index)
    SOB_dummies(:,n)=double(pob==SOB_index(n));
end
QTR1_3 = zeros(obs,3);
for q = 1:3
       QTR1_3(:,q)=double(qob==(q+1));
end
YOB_index = [31 32 33 34 35 36 37 38 39];
YR21_29 = zeros(obs,length(YOB_index));
for yy=1:length(YOB_index)
    YR21_29(:,yy)=double(yob==YOB_index(yy));
end
SOB_QTR_dummies=zeros(obs,length(SOB_index)*3);
for n=1:length(SOB_index)
    for q = 1:3
        idx = (n-1)*3 + q; %loop over quarters
        SOB_QTR_dummies(:,idx)=SOB_dummies(:,n).*QTR1_3(:,q);
    end
end

% QTR_YOB_POB_dummies
SOB_QTR_YOB_dummies=zeros(size(QTR1_3,1),size(YR21_29,2)*length(SOB_index)*size(QTR1_3,2));

for yy = 1:size(YR21_29,2)
    for n=1:length(SOB_index)
        for q = 1:size(QTR1_3,2)
            idx = (yy-1)*length(SOB_index)*size(QTR1_3,2)...
                + (n-1)*3 + q; %loop over years, states, quarters

            SOB_QTR_YOB_dummies(:,idx)=YR21_29(:,yy).*SOB_dummies(:,n).*QTR1_3(:,q);
        end
    end
end
SOB_YOB_dummies=zeros(size(QTR1_3,1),size(YR21_29,2)*length(SOB_index));
for yy = 1:size(YR21_29,2)
    for n=1:length(SOB_index)
            idx = (yy-1)*length(SOB_index)...
                + n; %loop over years, states

            SOB_YOB_dummies(:,idx)=YR21_29(:,yy).*SOB_dummies(:,n);
    end
end
QTR_YOB_dummies=zeros(size(QTR1_3,1),size(YR21_29,2)*3);
for yy = 1:size(YR21_29,2)
   for q = 1:size(QTR1_3,2)
   
            idx = (yy-1)*3 + q; %loop over years, qtr

            QTR_YOB_dummies(:,idx)=YR21_29(:,yy).*QTR1_3(:,q);
    end
end

% define (1530) instruments
%Z=[QTR1_3 QTR_YOB_dummies SOB_QTR_dummies SOB_YOB_dummies SOB_QTR_YOB_dummies];
%Z=[QTR1_3 QTR_YOB_dummies SOB_QTR_dummies SOB_YOB_dummies];
Z=[QTR1_3 QTR_YOB_dummies SOB_QTR_dummies];
%Z=[QTR_YOB_dummies SOB_QTR_dummies];
%Z=[QTR1_3 QTR_YOB_dummies];
%Z=[QTR1_3];
W=[ones(obs,1)];

% remove small groups
temp = sum(Z,1); zero_idx = find(temp >= 5); % remove multicollinearities and small groups
Z = Z(:,zero_idx);
temp = sum(W,1); zero_idx = find(temp >= 5); % remove multicollinearities
W = W(:,zero_idx);
m=size(Z,2); n=size(Z,1); disp(sprintf(' K and N are %d and %d',m,n))

%create projection matrix and helper matrices that depend on Z only
WW=W'*W; WZ=W'*Z; 
Z1=Z; Z=Z-W*WW^-1*WZ; % partial out constant from Z
Pw=W*WW^-1*W'; Wcontrol = Pw*control; % partial out constant from controls
ZZ=Z'*Z; ZZ_inv=ZZ^-1;
ZZ_inv2 = (ZZ_inv)^0.5;ZZZ_inv = Z*ZZ_inv; %n x k matrix
ZZZ_inv2 = Z*ZZ_inv2; %n x k matrix
H = sum(ZZZ_inv2.^2,2); %n x 1 vector of P_ii 
disp(sprintf(' Check balance: diagonal terms of P are btw %1.3f and %1.3f ', min(H),max(H)))
P = Z*ZZ_inv*Z'; % projection matrix
PP = P.*P;
PP_off = PP - diag(H.*H);
P_off = P - diag(H);
adj_cross = (1-H)*(1-H)' + PP; % cross fit variance need a different adjustment
PP_off_adj = PP_off./adj_cross;

kappa = sum(H.^2)/m; % for HLML heteroscedasticity-robust standard error
tau = m/n;

%theoretical first stage and F tilde (mu^2/\sqrt{K \Upsilon})
ZXmean = Z'*educbar;
pi = ZZ_inv*ZXmean;
Xmean = Z*pi;  

% f-statistic (linear case)
Xmean_l_demean = Xmean - mean(Xmean);
Xmean_var = ones(length(Xmean_l_demean),1);
Upsilon = 2*Xmean_var'*PP_off*Xmean_var/m;
Ftheo = Xmean_l_demean'*P_off*Xmean_l_demean/sqrt(m)/sqrt(Upsilon); 

%% Simulation block

sims = 1000

% Initialise variables:
F_SY_sims = zeros(sims,1); Fhat_sims = zeros(sims,1);
lAR_stat_sims = zeros(sims,1); lAR_sigma_sims = zeros(sims,1);
mKLM=zeros(grid_points,1);
%thli=zeros(grid_points,1);
%CLR=zeros(grid_points,1);
MCLR=zeros(grid_points,1);
AR_rej=zeros(grid_points,1);
% p_CLR=0; p_thli=0;
p_mKLM=0; p_mCLR=0; 

for i=1:grid_points
parfor s = 1:sims
    s
    rng(s)
    
    % simulate [Y,X] 
    v = normrnd(zeros(n,1),1); epsilon = normrnd(zeros(n,1),1);
    X = beta*Xmean + v;
    Y = beta*X + Wcontrol + v + beta*epsilon;

    % accumulate the outer product
    WX=W'*X;    WY=W'*Y;

    Y=Y-W*WW^-1*WY;    X=X-W*WW^-1*WX;
    Ybar = [Y,X];    YPY = Ybar'*P*Ybar;    Omega = Ybar'*Ybar - YPY;
    XY = X.*Y; YX = Y.*X; XX = X.*X; YY = Y.*Y;    XYYX = -XY-YX;
    YPY_off = Ybar'*P_off*Ybar; %leave-one-out version for JAK:

    % accumulate with P 
    ZX=Z'*X;   XZ=X'*Z;    ZY=Z'*Y;
    pi_hat = ZZ_inv*ZX;
    PX = ZZZ_inv*ZX;    MX = X - PX;
    MY = Y - ZZZ_inv*ZY;
    XMY = X.*MY; YMX = Y.*MX; XMX = X.*MX; YMY = Y.*MY;
    XMYYMX = -XMY-YMX; ZYbar = Z'*Ybar;

    % added for LR statistic computation 
    MYbar = Ybar - ZZZ_inv*ZYbar;
    H1=MYbar'*MYbar;
    
    %homoskedastic F hat
    v_hat = X - Z*pi_hat;
    F_SY = pi_hat'*ZZ*pi_hat/(v_hat'*v_hat/(n-m))/m;
    F_SY_sims(s,:) = F_SY;
   
    %% mKLM by Hansen_Hausman_Newey (2008, JBES)

    P2=(eye(n)-Pw)*Z*pinv(Z'*(eye(n)-Pw)*Z)*Z'*(eye(n)-Pw);
    gamma_hat=pinv(Wcontrol'*Wcontrol)*Wcontrol'*(Y-X*betagrid(i));
    M_full=eye(n)-Z*pinv(Z'*Z)*Z';
    H1_full = Ybar'*M_full*Ybar;

    b0=[1;-betagrid(i)];
    sigma_u0=b0'*H1_full*b0/(n-m-1);
    lam_tilde=((Y-X*betagrid(i))'*P2*(Y-X*betagrid(i)))/((Y-X*betagrid(i))'*M_full*(Y-X*betagrid(i)));
    Upsilon=P2*X; % (n*1)
    s_ee=(Y-X*betagrid(i))'*M_full*(Y-X*betagrid(i))/(n-m-1);
    s_ev=(Y-X*betagrid(i))'*M_full*X/(n-m-1);
    Y2_tilde=X-(Y-X*betagrid(i))*s_ev/s_ee;
    V_hat=M_full*Y2_tilde;
    kappa=sum(diag(P2).^2)/m;
    tau=m/n;
    Sigma_B=sigma_u0*(Y2_tilde'*P2*Y2_tilde+lam_tilde^2*Y2_tilde'*M_full*Y2_tilde);
    A_hat=(Upsilon'*(diag(P2)-tau))*(((Y-X*betagrid(i)-Wcontrol*gamma_hat).^2)'*V_hat/n);
    B_hat=m*(kappa-tau)*(V_hat'*diag(((Y-X*betagrid(i)-Wcontrol*gamma_hat).^2-sigma_u0))*V_hat)/(n*(1-2*tau+kappa*tau));
    Sigma_mKLM=Sigma_B+A_hat+A_hat'+B_hat;
    
    mKLM=(Y2_tilde'*P2*(Y-X*betagrid(i)))*pinv(Sigma_mKLM)*(Y2_tilde'*P2*(Y-X*betagrid(i)));

    if mKLM>mKLM95
        p_mKLM=p_mKLM+1;
    end

   
 %% Wald Test by Hausman et al (2012, QE)

    %Z_full = [W,Z];
    %M_full=eye(n)-Z_full*pinv(Z_full'*Z_full)*Z_full';
    %P_full = eye(n)-M_full;
    %P_h=P_full-diag(diag(P_full));
    %P_h=P-diag(diag(P));
    %M_h=eye(n)-P_h;
    
    %G_h=[Ybar,Wcontrol]'*P_h*[Ybar,Wcontrol];
    %H_h=[Ybar,Wcontrol]'*M_h*[Ybar,Wcontrol];
    
    %lam_h=min(eig(H_h^(-1/2)*G_h*H_h^(-1/2)));
    %alpha_h=min(eig(([Y,Z1]'*[Y,Z1])^(-1/2)*G_h*([Y,Z1]'*[Y,Z1])^(-1/2)));
    
    %X1=[X,Wcontrol];
    %beta_HLI=pinv(X1'*P_h*X1-lam_h*X1'*M_h*X1)*(X1'*P_h*Y-lam_h*X1'*M_h*Y);
    %beta_HLI=pinv(X'*P_h*X-alpha_h*(X'*X))*(X'*P_h*y-alpha_h*X'*y)
    %beta1_HLI=beta_HLI(1);
    
    %Q_h=X1'*(P_h-lam_h*M_h)*X1/n;
    %Q_h=X'*(P_h-alpha_h*eye(n))*X/n;
    %Q_inv_h=pinv(Q_h);
    
    %u_hat_h=Y-X1*beta_HLI;
    
    %Psi_HLI=Q_inv_h*...
    %(X1'*(P_h)*diag(u_hat_h.*u_hat_h)*(P_h)*X1+(X1.*u_hat_h)'*((P_h).^2)*(X1.*u_hat_h))/n...
    %*Q_inv_h;
    
    %t_hli = sqrt(n)*(beta1_HLI-betagrid(i))/sqrt((Psi_HLI(1,1)));

    %if t_hli>t_hli95
    %    p_thli=p_thli+1;
    %end

    %% CLR and M-CLR Statistics

    b0=[1;-betagrid(i)];
    A0=[betagrid(i);1];

    Omega_hat=H1/(n-m);
    lam=min(eig(H1^(-1/2)*YPY*H1^(-1/2)));

    S = ZZ_inv2*Z'*Ybar*b0*(b0'*Omega_hat*b0)^(-1/2);
    T = ZZ_inv2*Z'*Ybar*Omega_hat^(-1)*A0*(A0'*Omega_hat^(-1)*A0)^(-1/2);

    LR=S'*S-lam*(n-m);

    %c=c95(T,m);
    %if LR>c
    %p_CLR=p_CLR+1;
    %end

    c_m=c95_modified(T,n,m);
    if LR>c_m
    p_mCLR=p_mCLR+1;
    end

    %% Assemble Jaccknife AR statistic components:

    %covariance matrix terms:
    Sigma1_gg_0 = YMY'*PP_off_adj*YMY; 
    Sigma1_gg_1 = 2 * YMY'*PP_off_adj*XMYYMX; 
    Sigma1_gg_2 = 2 * XMX'*PP_off_adj'*YMY + XMYYMX'*PP_off_adj*XMYYMX; 
    Sigma1_gg_3 = 2 * XMX'*PP_off_adj*XMYYMX; 
    Sigma1_gg_4 = XMX'*PP_off_adj*XMX;

    % identification strength
    Fhat_sims(s,:) = (X'*P_off*X )/sqrt(m)/sqrt(2*Sigma1_gg_4/m);
 
    % assemble the statistic (functions)
     lAR_func = @(b) YPY_off(1,1) - 2 * YPY_off(1,2) * b + YPY_off(2,2)*b^2;
     Sigma1_gg_func = @(b)  2/m*(Sigma1_gg_0 + Sigma1_gg_1*b + Sigma1_gg_2*b^2 +...
        Sigma1_gg_3*b^3 + Sigma1_gg_4*b^4);

    % evaluate jackknife AR under the null
    lAR_stat_sims(s,:) = 1/sqrt(m)*lAR_func(betagrid(i));
    lAR_sigma_sims(s,:) = Sigma1_gg_func(betagrid(i));
end
    AR_rej1 = (lAR_stat_sims./sqrt(lAR_sigma_sims)) > JAK95;
    mKLM(i)=p_mKLM/sims;
    MCLR(i)=p_mCLR/sims;
    AR_rej(i) = mean(AR_rej1);
    p_mKLM=0; p_mCLR=0; 
end

%% Save figure
 
figurename = strcat('power_0.25percent_final', '.eps');
fig = figure(1);
betagrid_plot = betagrid;
plot(-betagrid_plot',mKLM(ismember(betagrid,betagrid_plot)),'-.',...
    -betagrid_plot,MCLR(ismember(betagrid,betagrid_plot))', '-.',...
    -betagrid_plot,AR_rej(ismember(betagrid,betagrid_plot))', '-.',...
    'LineWidth',2);

legend({sprintf('mKLM'),sprintf('MCLR'),sprintf('Jackknife AR')},...
    'Location','southoutside','Orientation','horizontal','FontSize',14 );
legend('boxoff');
xlabel('\beta'); xlim([min(betagrid_plot) max(betagrid_plot)])
ylabel('Rejection Probability of 5% Test, H_0:\beta=0.1','FontSize',14 );
ylim([0 1])
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,figurename,'-depsc','-r0')

