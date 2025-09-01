% here i am trying to check if the F-statistic can be used to evaluate
% significance of ISC on iid training data. 
clear all
T=20; N=2;
df3=N*(T-1);
df2=T*(N-1); 
df1=T-1; 
Nrand=100000;

for i=Nrand:-1:1

    x = randn(T,N);
    
    % equalizing means makes both measures the same
    % x = x - repmat(mean(x),T,1);

    rw(i) = sum(var(x))*(T-1);
    rt(i) = N^2*var(mean(x,2))*(T-1);
    rb(i) = rt(i) - rw(i);
    rho = 1/(N-1)*rb(i)/rw(i); % ISC
    
    Scca(i) = (rho+1/(N-1))./(1-rho); 
    
    Fcca(i) = (N-1)*rho+1; % incorrectly assumes that rt and rw are independent
    
    sw = sum(var(x,[],2))*(N-1);
    st = var(x(:))*(T*N-1);
    sb = st - sw;
    Slda(i) = sb/sw;    


end

subplot(2,2,1); 
plot(Scca,Slda,'.')
xlabel('S estimated from \rho')
ylabel('S = s_b/s_w')

subplot(3,2,2); 
[Prand,fbin] = hist(Fcca,100); 
P = fpdf(fbin,df1,df3);
bar(fbin,Prand/Nrand/mean(diff(fbin))); hold on
plot(fbin,P,'r'); hold off
xlabel('new F computed from \rho')

subplot(3,2,4); 
F = df2/df1*Scca; % only correct if means as equalized 
[Prand,fbin] = hist(F,100); 
P = fpdf(fbin,df1,df2);
bar(fbin,Prand/Nrand/mean(diff(fbin))); hold on
plot(fbin,P,'r'); hold off
xlabel('F computed from \rho')

subplot(4,2,5); 
[Prand,fbin] = hist(rw,100); 
P = chi2pdf(fbin,df3);
bar(fbin,Prand/Nrand/mean(diff(fbin))); hold on
plot(fbin,P,'r'); hold off
xlabel('r_W')

subplot(4,2,7);
[Prand,fbin] = hist(rt/N,100); 
P = chi2pdf(fbin,df1);
bar(fbin,Prand/Nrand/mean(diff(fbin))); hold on
plot(fbin,P,'r'); hold off
xlabel('r_T')


subplot(3,2,6); 
F = df2/df1*Slda;
[Prand,fbin] = hist(F,100); 
P = fpdf(fbin,df1,df2);
bar(fbin,Prand/Nrand/mean(diff(fbin))); hold on
plot(fbin,P,'r'); hold off
xlabel('F computed from S')
