function c_95m=c95_modified(T,n,k)

TT=T'*T;

M=5000;
LRm=zeros(M,1);

for i=1:M

    
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

LRm(i,1)=(n-k)*(SS/S_star2-lam_hat);
end

LRm=sort(LRm);
c_95m=LRm(length(LRm)*0.95);