function [t0,t] = Opt_P1_k(gam,k)



A=sum(gam);
syms z
eqn=z*log(z)-z+1==A;
optz=double(solve(eqn,z));
t0=(optz-1)/(A+optz-1);

for i=1:k


t(i)=gam(i)/(A+optz-1);

end
