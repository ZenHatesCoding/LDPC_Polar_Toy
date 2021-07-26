clc;
clear;
close all; 

crcg = fliplr([1 1 1 0 0 0 1 0 0 0 0 1]);
crcL = length(crcg)-1;

List_decoding_or_Not = 1;

N = 128;
load('Polar_reliability.mat');
Q1 = Q(Q<=N); % reliability sequence for N
% Q1 = [1 2 3 5 9 4 6 10 7 11 13 8 12 14 15 16];

K = 119; 
if List_decoding_or_Not
    A = K - crcL;
end
Nblock = 100;


F = Q1(1:(N-K)); % frozen position Q1(1:(N-K))
                 % Message position Q1(N-K+1:end)

if List_decoding_or_Not
    msg = randi([0 1],1,A*Nblock);
else
    msg = randi([0 1],1,K*Nblock);
end
% Encoding
u = [];
for idx = 1:Nblock
    if List_decoding_or_Not
        u_block = Polar_Encoder_wCRC(msg((idx-1)*A+1:idx*A),N,A,Q1,crcg);
    else
        u_block = Polar_Encoder(msg((idx-1)*K+1:idx*K),N,K,Q1);
    end
    u = [u, u_block];
end

% transmission
EbNo_dB = 6;
if List_decoding_or_Not
    Rate = A/N;
else
    Rate = K/N;
end
EbNo = 10^(EbNo_dB/10);
sigma = sqrt(1/(2*Rate*EbNo));


s = 1-2*u;
r = s + sigma*randn(1,length(s));
msgcap = [];
NL = 4;
for idx = 1:Nblock
    if List_decoding_or_Not
        msgcap_block = SSCList_Polar_Decoder(r((idx-1)*N+1:idx*N),N,Q1,F,NL,crcg);
    else
        msgcap_block = SSC_Polar_Decoder(r((idx-1)*N+1:idx*N),N,Q1,F);
    end
    msgcap = [msgcap msgcap_block];
end
disp(['BER w/o polar code:',num2str((sum(abs(u-(r<0)))/length(msg)),'%e')]);
disp(['BER w polar code:',num2str((sum(abs(msg-msgcap))/length(msg)),'%e')]);

