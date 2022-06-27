clear all
clc
m=0;




a=[2 2.5 3 3.5 4];

gam1=10^(2.2);
gam2=10.^([10 7 4 1 -2]./10)';
T=[ 2.6810    3.5954    4.9119    6.7964    9.5953];
  


syms z
for i=1:length(a)
A=gam1+gam2(i);

eqn=z*log(z)-z+1==A;
optz(i,1)=double(solve(eqn,z));

 end   
opt0=(optz-1)./(A+optz-1);
opt1=gam1./(A+optz-1);
opt2=gam2./(A+optz-1);
R1_op=opt1.*log2(1+gam1.*(opt0./opt1));
R2_op=opt2.*log2(1+gam2.*(opt0./opt2));
R=R2_op./R1_op; 

hold on
semilogy(a,T,'o-')
semilogy(a,R,'*-')

grid on
hold off