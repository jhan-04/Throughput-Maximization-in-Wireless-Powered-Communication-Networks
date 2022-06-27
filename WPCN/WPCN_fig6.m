clear all
clc

gam1=10^(2.2);
gam2=10.^([10 7 4 1 -2]./10)';


a=[2 2.5 3 3.5 4]';
n=1000;
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


semilogy(a,R,'*-')

grid on


%TDMA

a=[2 2.5 3 3.5 4]';
gam1_t=10^(1.3);
gam2_t=10.^([7.1 5.5 4 2.4 0.9]./10)';
for i=1:length(a)
A=gam1_t+gam2_t(i);
eqn=z*log(z)-z+1==A;
optz_t(i,1)=double(solve(eqn,z));
 end   
    
  
opt0_t=(optz_t-1)./(A+optz_t-1);
opt1_t=gam1_t./(A+optz_t-1);
opt2_t=gam2_t./(A+optz_t-1);
R1_op_t=opt1_t.*log2(1+gam1_t.*(opt0_t./opt1_t));
R2_op_t=opt2_t.*log2(1+gam2_t.*(opt0_t./opt2_t));
R_t=R2_op_t./R1_op_t; 
hold on

semilogy(a,R_t,'-o')

xlabel('Pathloss Exponent a')
ylabel('R2(tou*)/R1(tou*)')

legend('WPCN','Conventional TDMA network')
hold off