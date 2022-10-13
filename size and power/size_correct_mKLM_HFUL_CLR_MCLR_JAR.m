function [mKLM95, t_hli95, CLR095, MCLR095, JAK95]=size_correct_mKLM_HFUL_CLR_MCLR_JAR(n,m,k,rho,delta2,beta10)

% Initialise vectors:
mKLM=zeros(m,1);
t_hli=zeros(m,1);
CLR=zeros(m,1);
MCLR=zeros(m,1);
JAK=zeros(m,1);

% True beta:
%beta10=0;
beta20=0;
%beta0=[beta10;beta20];

for j=1:m

    % Create instruments, and fix them:
    k1=1;
    k2=k-k1;
    %Z1=randn(n,k1);
    Z1=ones(n,k1);
    %Z2=randn(n,k2);
    
    Z20=randn(n,1);
    %Z2=[Z20, Z20.^2, Z20.^3, Z20.^4, randn(n,k2-4)];
    Z2=[Z20, Z20.^2, Z20.^3, randn(n,k2-3)];
    %Z2=[Z20, Z20.^2, randn(n,k2-2)];
    
    Z=[Z1,Z2];
        
    % Error variables:
    u=randn(n,1); % Homoskedasticity
    v2=(1-rho^2)^0.5*randn(n,1)+rho*u/(1/1)^0.5; %construct random error terms, which are a function of u (degree of endogeneity)
    
    % PUT THE MODEL TOGETHER
    pi2=(delta2/sum(sum(Z2'*Z2-Z2'*Z1*pinv(Z1'*Z1)*Z1'*Z2)))^0.5; % de-meaned coefficient, determined by delta^2 (order of magnitude)
    Y2=Z*(repmat(pi2,k,1))+v2; %endogenous regressor, since it has the v2 element in it
    y1=beta20*Z1+beta10*Y2+u; %reduced form of y1 with intercept
    Y=[y1,Y2]; % (n x 2) matrix, when l=1 
    
    % Transformations of instrument matrices
    ZZ=Z'*Z; ZZ_inv=ZZ^-1;
    ZZ_inv2 = (ZZ_inv)^0.5;
    ZZZ_inv = Z*ZZ_inv; %n x k matrix
    ZZZ_inv2 = Z*ZZ_inv2; %n x k matrix
    P_diag = sum(ZZZ_inv2.^2,2); %n x 1 vector of P_ii 

    % Projection matrices + related transformations
    P=Z*pinv(Z'*Z)*Z';
    PP = P.*P;
    PP_off = PP - diag(P_diag.*P_diag);
    P_off = P - diag(P_diag);
    M=eye(n)-P;
    P1=Z1*pinv(Z1'*Z1)*Z1';
    P2=(eye(n)-P1)*Z2*pinv(Z2'*(eye(n)-P1)*Z2)*Z2'*(eye(n)-P1);
    G=Y'*(P-P1)*Y;
    H=Y'*M*Y;

    % ------------------------------------------------------------------------
    % mKLM by Hansen_Hausman_Newey (2008, JBES)

    gamma_hat=pinv(Z1'*Z1)*Z1'*(y1-Y2*beta10);

    b0=[1;-beta10];
    sigma_u0=b0'*H*b0/(n-k);
    lam_tilde=((y1-Y2*beta10)'*P2*(y1-Y2*beta10))/((y1-Y2*beta10)'*M*(y1-Y2*beta10));
    Upsilon=P2*Y2; % (n*1)
    s_ee=(y1-Y2*beta10)'*M*(y1-Y2*beta10)/(n-k);
    s_ev=(y1-Y2*beta10)'*M*Y2/(n-k);
    Y2_tilde=Y2-(y1-Y2*beta10)*s_ev/s_ee;
    V_hat=M*Y2_tilde;
    kappa=sum(diag(P2).^2)/k2;
    tau=k2/n;
    Sigma_B=sigma_u0*(Y2_tilde'*P2*Y2_tilde+lam_tilde^2*Y2_tilde'*M*Y2_tilde);
    A_hat=(Upsilon'*(diag(P2)-tau))*(((y1-Y2*beta10-Z1*gamma_hat).^2)'*V_hat/n);
    B_hat=k2*(kappa-tau)*(V_hat'*diag(((y1-Y2*beta10-Z1*gamma_hat).^2-sigma_u0))*V_hat)/(n*(1-2*tau+kappa*tau));
    Sigma_mKLM=Sigma_B+A_hat+A_hat'+B_hat;
    
    mKLM(j,1)=(Y2_tilde'*P2*(y1-Y2*beta10))*pinv(Sigma_mKLM)*(Y2_tilde'*P2*(y1-Y2*beta10));
    
    % ------------------------------------------------------------------------
    % Wald Test by Hausman et al (2012, QE)

    P_h=P-diag(diag(P));
    M_h=eye(n)-P_h;
    
    G_h=[Y,Z1]'*P_h*[Y,Z1];
    H_h=[Y,Z1]'*M_h*[Y,Z1];
    
    lam_h=min(eig(H_h^(-1/2)*G_h*H_h^(-1/2)));
    %alpha_h=min(eig(([Y,Z1]'*[Y,Z1])^(-1/2)*G_h*([Y,Z1]'*[Y,Z1])^(-1/2)));
    
    X=[Y2,Z1];
    beta_HLI=pinv(X'*P_h*X-lam_h*X'*M_h*X)*(X'*P_h*y1-lam_h*X'*M_h*y1);
    %beta_HLI=pinv(X'*P_h*X-alpha_h*(X'*X))*(X'*P_h*y-alpha_h*X'*y)
    beta1_HLI=beta_HLI(1);
    
    Q_h=X'*(P_h-lam_h*M_h)*X/n;
    %Q_h=X'*(P_h-alpha_h*eye(n))*X/n;
    Q_inv_h=pinv(Q_h);
    
    u_hat_h=y1-X*beta_HLI;
    
    Psi_HLI=Q_inv_h*...
    (X'*(P_h)*diag(u_hat_h.*u_hat_h)*(P_h)*X+(X.*u_hat_h)'*((P_h).^2)*(X.*u_hat_h))/n...
    *Q_inv_h;
    
    t_hli(j,1) = sqrt(n)*(beta1_HLI-beta10)/sqrt((Psi_HLI(1,1)));


    %prod_mat = (Y'*Y)^(-1)*(Y'*P_off*Y);
    %alpha_tilde = min(eig(prod_mat));
    %alpha_hat = (alpha_tilde - (1-alpha_tilde)/n)/(1-(1-alpha_tilde)/n);
    %H_hat = Y2'*P_off*Y2-alpha_hat*Y2'*Y2;
    %beta_liml_hat = (H_hat)^(-1)*(Y2'*P_off*y1-alpha_hat*Y2'*y1);

    %epsilon_hat = y1 - beta_liml_hat*Y2;
    %gamma_hat = (Y2'*epsilon_hat)/(epsilon_hat'*epsilon_hat);
    %X_hat = Y2 - epsilon_hat*gamma_hat';
    %X_dot = P*X_hat;
    %sigma_hat1 = epsilon_hat'*(X_dot'*X_dot - X_dot'*diag(P_diag)*X_hat - X_hat'*diag(P_diag)*X_dot)*epsilon_hat;
    %sigma_hat2 = (epsilon_hat'*(Z*Z'*X_hat))'*(epsilon_hat'*(ZZZ_inv*ZZZ_inv'*X_hat));
    %sigma_hat = sigma_hat1 + sigma_hat2;
    %var_sigma_hat = H_hat^(-1)*sigma_hat*H_hat^(-1);

    %t_hli(j,1) = (beta_liml_hat-beta10)/sqrt(var_sigma_hat);

    % ------------------------------------------------------------------------
    % CLR by Moreira (2003)
    
    A0=[beta10;1];

    % Estimate of Omega
    Omega_hat=H/(n-k);
    
    % Define T_bar
    T=ZZ_inv2*Z'*Y*Omega_hat^(-1)*A0*(A0'*Omega_hat^(-1)*A0)^(-1/2);
    TT=T'*T;
    
    Q_1=chi2rnd(1);
    Q_2=chi2rnd(k-1);
    CLR(j,1)=0.5*(Q_1+Q_2-TT+sqrt((Q_1+Q_2+TT)^2-4*Q_2*(TT)));
    
    % ------------------------------------------------------------------------
    % Modified CLR
        
    % Sufficient statistics:
    mu = zeros(k,1);
    Ik = eye(k);
    S=(mvnrnd(mu,Ik))';
    W=wishrnd(eye(2),n-k);
    
    SS=S'*S;
    ST=S'*T;
    S_star2 = W(1,1);
    
    GG=[SS,ST;ST',TT];
    
    lam_hat=min(eig(W^(-1/2)*GG*W^(-1/2)));
    
    MCLR(j,1)=(n-k)*(SS/S_star2-lam_hat);   
    % ------------------------------------------------------------------------
    % Mikusheva and Sun (2021) Jackknife AR
    
    % Cross fit variance needs an adjustment relative to standard AR
    adj_cross = (1-P_diag)*(1-P_diag)' + PP;
    PP_off_adj = PP_off./adj_cross;

    YPY = Y'*P_off*Y;
    YPY1 = YPY(1,1); YPY2 = YPY(1,2); YPY3 = YPY(2,2);
    
    ZX=Z'*Y2; ZY=Z'*y1; 
    PX = ZZZ_inv*ZX;    MX = Y2 - PX;
    MY = y1 - ZZZ_inv*ZY;
    XMY = Y2.*MY; YMX = y1.*MX; XMX = Y2.*MX; YMY = y1.*MY;
    XMYYMX = -XMY-YMX;
    
    Sigma_gg_0 = YMY'*PP_off_adj*YMY; 
    Sigma_gg_1 = 2 * YMY'*PP_off_adj*XMYYMX; 
    Sigma_gg_2 = 2 * XMX'*PP_off_adj'*YMY + XMYYMX'*PP_off_adj*XMYYMX; 
    Sigma_gg_3 = 2 * XMX'*PP_off_adj*XMYYMX; 
    Sigma_gg_4 = XMX'*PP_off_adj*XMX;
    
    lAR_func = @(b) 1/sqrt(k)*(YPY1 - 2 * YPY2 * b + YPY3*b.^2); 
    Sigma_gg_func = @(b)  2/k*(Sigma_gg_0 + Sigma_gg_1*b + Sigma_gg_2*b.^2 +...
                Sigma_gg_3*b.^3 + Sigma_gg_4*b.^4);
    JAK(j,1) = (lAR_func(beta10)./sqrt(Sigma_gg_func(beta10)));

end

% ------------------------------------------------------------------------

mKLM=sort(mKLM);
t_hli = sort(t_hli);
CLR=sort(CLR);
MCLR=sort(MCLR);
JAK=sort(JAK);

mKLM95=mKLM(m*0.95);
t_hli95=t_hli(m*0.95);
CLR095=CLR(m*0.95);
MCLR095=MCLR(m*0.95);
JAK95=JAK(m*0.95);
