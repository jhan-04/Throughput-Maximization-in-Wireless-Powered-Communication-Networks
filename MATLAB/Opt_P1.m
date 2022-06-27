function [t] = Opt_P1(gam)
gam1=gam(1);
gam2=gam(2);
A=gam1+gam2;
syms z
eqn=z*log(z)-z+1==A;
optz=double(solve(eqn,z));
opt0=(optz-1)/(A+optz-1);
opt1=gam1/(A+optz-1);
opt2=gam2/(A+optz-1);
R1_op=opt1.*log2(1+gam1*(opt0./opt1));
R2_op=opt2.*log2(1+gam2*(opt0./opt2));
R_op=R1_op+R2_op;
t=[opt0 opt1 opt2];
end
