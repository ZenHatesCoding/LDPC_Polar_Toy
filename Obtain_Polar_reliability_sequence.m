TxT = fopen('Polar_reliability_sequence_raw.txt');
A = fscanf(TxT,'%d',[16 128]).';
fclose(TxT);

A2 = A(:,2:2:end);
Q = A2(:)+1;

save('Polar_reliability.mat','Q');