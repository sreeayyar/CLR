function [p_CLR, p_mCLR]=criticalval(n,m,k,tau)

p_CLR=0;
p_mCLR=0;

T=(sqrt(tau/k)*ones(k,1));

for j=1:m

    %CLR
    c=c95(T,k);
    p_CLR=p_CLR+c;
    
    %modified CLR
    c_m=c95_modified(T,n,k);
    p_mCLR=p_mCLR+c_m;

end

% evaluate the empirical sizes
p_CLR=p_CLR/m;
p_mCLR=p_mCLR/m;
