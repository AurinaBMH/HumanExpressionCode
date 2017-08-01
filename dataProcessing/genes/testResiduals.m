% test if residuals look fine
figure;
for j=1:10
y = randsample(1626900,1000, true);
A = distExpVect(y,1); B = Residuals(y,1); 
%Residuals(y);
subplot(2,5,j); plot(A, B,'.k');
[r,p] = corr(A,B)
end

