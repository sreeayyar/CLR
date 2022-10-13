function c95=c95(T,k)

TT=T'*T;

M=5000;
LR=zeros(M,1);

for i=1:M
Q_1=chi2rnd(1);
Q_2=chi2rnd(k-1);
LR(i,1)=0.5*(Q_1+Q_2-TT+sqrt((Q_1+Q_2+TT)^2-4*Q_2*(TT)));
end
LR=sort(LR);
c95=LR(length(LR)*0.95);

