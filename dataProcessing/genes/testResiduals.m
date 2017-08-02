% test if residuals look fine
figure;
for j=1:10
y = randsample(1626900,1000, true);
A = distExpVect(y,1); B = Residuals(y,1); 
%Residuals(y);
[r,p] = corr(A,B);
subplot(2,5,j); plot(A, B,'.k'); title(sprintf('r=%d, p=%d',r,p)); 

end

