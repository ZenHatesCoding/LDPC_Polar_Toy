close all;
clear
num_block_K = 119;
m = 100;
Xin = randi([0,1],m,num_block_K);

if num_block_K == 119
    Hamming_matrix;
else
    [Hamming_H,Hamming_G] = hammgen(7);
end

num_block_M = size(Hamming_G,2);

Cin = mod(Xin*Hamming_G,2);

Cin = reshape(Cin.',1,[]);
x_Tx = 1-2*Cin; % [0,1] -> [1, -1]

x_Rx = x_Tx;

% error_idx = [ 6 103];
% x_Rx(error_idx) = -0.2*x_Rx(error_idx);
SNR = 8;
noise_power = 10*log10((1/(10^(SNR/10))));
x_Rx = x_Rx + wgn(1,length(x_Rx),noise_power,'dBW','real');


% SISO decision decoding
disp('SD Decoding:');
figure;
stem(x_Tx(1:num_block_M)); hold on; stem(x_Rx(1:num_block_M));

nlayer = 1;
niter = 16;
Offset = 1e-2;

x_Rx = rescale(x_Rx,-1,1);
C = [];
for idx = 1:m
%     C_temp = SISO_minSum_layered_decoder(x_Rx((idx-1)*num_block_M+1:idx*num_block_M),...
%         Hamming_H,nlayer,niter,Offset);
%     C = [C, C_temp<0];
    C_temp = Hamming_SIHO_Chase_Decoder(x_Rx((idx-1)*num_block_M+1:idx*num_block_M),...
             Hamming_H,3,10);
    C = [C, C_temp];
end
nEr1 = sum(abs(Cin-(x_Rx<0)));
nEr2 = sum(abs(Cin-(C)));
disp(['BER w/o FEC:',num2str(nEr1/length(Cin),'%e')]);
disp(['BER with FEC:',num2str(nEr2/length(Cin),'%e')]);

% Hard decision decoding
disp('HD Decoding:');
nEr1 = sum(abs(Cin-(x_Rx<0)));
x_Rx_new = [];
for idx = 1:m
    x_Rx_block = x_Rx((idx-1)*num_block_M+1:idx*num_block_M);
    x_Rx_block = Hamming_syndrome_decoder(x_Rx_block,Hamming_H);
%     syndrome = mod(Hamming_H*(x_Rx_block<0).',2);
%     [~,I] = ismember(Hamming_H.',syndrome.','rows');
%     jdx = find(I==1);
%     x_Rx_block(jdx) = -x_Rx_block(jdx);
    x_Rx_new = [x_Rx_new, x_Rx_block];
end
nEr2 = sum(abs(Cin-x_Rx_new));
disp(['BER w/o FEC:',num2str(nEr1/length(Cin),'%e')]);
disp(['BER with FEC:',num2str(nEr2/length(Cin),'%e')]);

