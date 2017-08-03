% test noise levels
fileNoise = 'PACall.csv';
noise = csvread(fileNoise);

[~,probeList] = intersect(noise(:,1),ProbeID, 'stable');
A = noise(probeList,2:end);
B = sum(A,2);

figure; histogram(B/size(A,2), 20); title('Histogram of signal in probes'); xlabel('proportion of signal in the probe'); ylabel('count');
figure; imagesc(A); title('Origianl binary signal representation'); xlabel('Sample'); ylabel('probe');


[~,~,probeList2] = intersect(noise(:,1),ProbeID, 'stable');
B = sum(A,2)/size(A,2); C = B>0.5; D = probeList2(C==1); E = A(D,:);
figure; imagesc(E);  title('Filtered binary signal representation'); xlabel('Sample'); ylabel('probe');


