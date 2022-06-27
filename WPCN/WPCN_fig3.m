%%when K is 0,ploting throughput ,Fig 3.
i=1;
gam1=10;
t0=[0:0.01:1];
t1=1-t0;
R1=t1.*log2(1+gam1*(t0./t1));
axis([0,1,0,2])
hold on
plot(t0,R1)
[M,I] = max(R1(:));
t0(I)
hold on
line([t0(I) t0(I)],ylim)
xlabel('t0')
ylabel('Throughput(bps/Hz)')
