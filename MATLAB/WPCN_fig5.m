clear all
clc
m=0;
mt1=0;
mt2=0;

gam1=10^(2.2);
gam2=10^(1);

    k=1000;
axis([0,1,0,1])
for c=1:1:k
    for d=1:1:k-c
        t1(c,d)=d/k;
        t2(c,d)=c/k;
        t0(c,d)=1-d/k-c/k;
     end
end    
    R1=t1.*log2(1+gam1*(t0./t1));
    R2=t2.*log2(1+gam2*(t0./t2));
    R=R1+R2;
    hold on

contour(t1,t2,R,8)
max(max(R))
%optimal
A=gam1+gam2;
syms z
eqn=z*log(z)-z+1==A;
optz=double(solve(eqn,z));
opt0=(optz-1)/(A+optz-1);
opt1=gam1/(A+optz-1);
opt2=gam2/(A+optz-1);
R1_op=opt1.*log2(1+gam1*(opt0./opt1));
R2_op=opt2.*log2(1+gam2*(opt0./opt2));
R_op=R1_op+R2_op
    
    plot(opt1,opt2,'*')

hold off