% comments
% polar code not working for low SNR (coding results in even worse performance)
% but for high SNR it achieves error free transmission
%%
close all;
clear

m = 256;
nC = 2; % number of channel
SE = 2;

LBC_or_Polar = 1; % 1 LBC (LN block code) 2 & polar code

Enable_WL_interleaving = 1;
num_block_interleave = 16;
switch LBC_or_Polar
    case 1
        Enable_SD_decoder = 1;
            SD_decoder_type = 2; % 1 minSum decoder 2 chase II decoder
        Enable_HD_decoder = 1;
        num_block_K = 119;

        Xin = randi([0,1],m*SE*nC,num_block_K);
    case 2
        Enable_List_decoding = 1;
        NL = 4; % used in list decoder only
        crcg = fliplr([1 1 1 0 0 0 1 0 0 0 0 1]);
        crcL = length(crcg)-1;

        N = 128;
        load('Polar_reliability.mat');
        Q1 = Q(Q<=N); % reliability sequence for N
        
        K = 119; 
        if Enable_List_decoding
            A = K - crcL;
        end

        Nblock = m*nC*SE;
        F = Q1(1:(N-K)); % frozen position Q1(1:(N-K))
                         % Message position Q1(N-K+1:end)
end

%% Encoding
switch LBC_or_Polar
    case 1
        Hamming_matrix;
        num_block_M = size(Hamming_G,2);

        Cin = mod(Xin*Hamming_G,2);

        Cin = reshape(Cin.',1,[]);
        
        
        Cin_mat = reshape(Cin,[],nC).'; 

    case 2  
        if Enable_List_decoding
            msg = randi([0 1],1,A*Nblock);
        else
            msg = randi([0 1],1,K*Nblock);
        end
        u = [];
        for idx = 1:Nblock
            if Enable_List_decoding
                u_block = Polar_Encoder_wCRC(msg((idx-1)*A+1:idx*A),N,A,Q1,crcg);
            else
                u_block = Polar_Encoder(msg((idx-1)*K+1:idx*K),N,K,Q1);
            end
            u = [u, u_block];
        end
        Cin = u;
        Cin_mat = reshape(Cin,[],nC).'; 
        msg_mat = reshape(msg,[],nC).';
end
%% PAM4 modulation 

x_Tx = tx_PAM4_mod(~Cin);


%% channel noise loading
SNR_good = 14;
SNR = [SNR_good 14 14 14];
SNR = SNR(1:nC);
noise_power = 10.*log10((1./(10.^(SNR./10))));
x_Tx_mat = reshape(x_Tx,num_block_interleave*nC,[]);
if Enable_WL_interleaving
    x_Tx_wnoise_wdm = zeros(nC,length(x_Tx)/nC);    
    for idx = 1:nC
        x_Tx_wnoise_wdm(idx,:) = reshape(x_Tx_mat((idx-1)*num_block_interleave+1:idx*num_block_interleave,:),1,[]);
    end
else
    x_Tx_wnoise_wdm = reshape(x_Tx,[],nC).';
end
for idx = 1:nC
    x_Tx_wnoise_wdm(idx,:) = x_Tx_wnoise_wdm(idx,:)...
                             + wgn(1,size(x_Tx_wnoise_wdm,2),...
                              noise_power(idx),'dBW','real');
 
end

%% decoding
switch LBC_or_Polar
    case 1
        if Enable_SD_decoder
            % % SISO decision decoding
            disp('SD Decoding:');

            nlayer = 3;
            niter = 1;
            Offset = 1e-2;
            L_Rx_wdm = zeros(size(x_Tx_wnoise_wdm,1),size(x_Tx_wnoise_wdm,2)*SE);
            
            for idx = 1:nC
                L_Rx_wdm(idx,:) = compute_LLR(x_Tx_wnoise_wdm(idx,:),...
                                  SE,10^(SNR(idx)/10));           
            end
            
            x_Rx_wdm = x_Tx_wnoise_wdm;
            
            if Enable_WL_interleaving
                L_Rx_mat = [];
                x_Rx_mat = [];
                
                for idx = 1:nC
                   L_Rx_mat = [L_Rx_mat;...
                               reshape(L_Rx_wdm(idx,:),num_block_interleave*SE,[])]; 
                   x_Rx_mat = [x_Rx_mat;...
                               reshape(x_Tx_wnoise_wdm(idx,:),num_block_interleave,[])];
                end
                Lblc_size_temp = size(L_Rx_mat,2)/nC;
                xblc_size_temp = size(x_Rx_mat,2)/nC;
                for idx = 1:nC
                    L_Rx_wdm(idx,:) = reshape(L_Rx_mat(:,(idx-1)*Lblc_size_temp+1:idx*Lblc_size_temp),1,[]);
                    x_Rx_wdm(idx,:) = reshape(x_Rx_mat(:,(idx-1)*xblc_size_temp+1:idx*xblc_size_temp),1,[]);
                end

            end
            % channels
            nEr1 = [];
            nEr2 = [];
%             L_Rx_wdm_afterSD = zeros(size(L_Rx_wdm));
            for jdx = 1:nC
                x_Rx = x_Rx_wdm(jdx,:);
                L_Rx = L_Rx_wdm(jdx,:);
                % trick
                L_Rx = L_Rx-mean(L_Rx);
                L_Rx = pwr_normalization(L_Rx);

                C = [];
                
                for idx = 1:m*SE
                    L_Rx_block = L_Rx((idx-1)*num_block_M+1:idx*num_block_M);
                    switch SD_decoder_type
                        case 1
                            C_temp = SISO_minSum_layered_decoder(L_Rx_block,...
                                Hamming_H,nlayer,niter,Offset);
                            C = [C, C_temp<0];
                        case 2
                            C_temp = Hamming_SIHO_Chase_Decoder(L_Rx_block,...
                                Hamming_H,3,10);
                            C = [C, C_temp];
                    end

                end
                C_rx = rx_PAM4_Decode(x_Rx);
                C_rx2 = ~C_rx;
                nEr1 = [nEr1,sum(abs(Cin_mat(jdx,:)-(C_rx2)))];
                nEr2 = [nEr2,sum(abs(Cin_mat(jdx,:)-(~C)))];
            end

            nEr1 = [nEr1, mean(nEr1)];
            nEr2 = [nEr2, mean(nEr2)];

            disp('BER w/o SD Hamming code:');
            nEr1/size(Cin_mat,2)
            disp('BER with SD Hamming code:');
            nEr2/size(Cin_mat,2)
        end

        if Enable_HD_decoder
            % Hard decision decoding
            disp('HD Decoding:');
            x_Rx_wdm = x_Tx_wnoise_wdm;
            
            if Enable_WL_interleaving           
                x_Rx_mat = [];               
                for idx = 1:nC
                    x_Rx_mat = [x_Rx_mat;...
                               reshape(x_Tx_wnoise_wdm(idx,:),num_block_interleave,[])];
                end
                xblc_size_temp = size(x_Rx_mat,2)/nC;
                for idx = 1:nC
                    x_Rx_wdm(idx,:) = reshape(x_Rx_mat(:,(idx-1)*xblc_size_temp+1:idx*xblc_size_temp),1,[]);
                end
            end
            
            
            % channel
            nEr1 = [];
            nEr2 = [];
            for jdx = 1:nC
                C_rx = rx_PAM4_Decode(x_Rx_wdm(jdx,:));
                C_rx2 = ~C_rx;
                nEr1 = [nEr1,sum(abs(Cin_mat(jdx,:)-(C_rx2)))];
                C_rx_BPSK = 1-2*C_rx2;
                x_Rx_new = [];
                
%                 C_rx_BPSK = L_Rx_wdm_afterSD(jdx,:)<0;
                for idx = 1:m*SE
                    x_Rx_block = C_rx_BPSK((idx-1)*num_block_M+1:idx*num_block_M);
                    % capable of correcting one error
                    x_Rx_block = Hamming_syndrome_decoder(x_Rx_block,Hamming_H);
                    x_Rx_new = [x_Rx_new, x_Rx_block];
                end
                nEr2 = [nEr2, sum(abs(Cin_mat(jdx,:)-x_Rx_new))];
            end
            nEr1 = [nEr1, mean(nEr1)];
            nEr2 = [nEr2, mean(nEr2)];

%             disp('BER w/o HD Hamming code:');
%             nEr1/size(Cin_mat,2)
            disp('BER with HD Hamming code:');
            nEr2/size(Cin_mat,2)
        end
    case 2
        L_Rx_wdm = zeros(size(x_Tx_wnoise_wdm,1),size(x_Tx_wnoise_wdm,2)*SE);
        for idx = 1:nC
            L_Rx_wdm(idx,:) = compute_LLR(x_Tx_wnoise_wdm(idx,:),...
                              SE,10^(SNR(idx)/10));           
        end

        x_Rx_wdm = x_Tx_wnoise_wdm;
        if Enable_WL_interleaving
            L_Rx_mat = [];
            x_Rx_mat = [];

            for idx = 1:nC
               L_Rx_mat = [L_Rx_mat;...
                           reshape(L_Rx_wdm(idx,:),num_block_interleave*SE,[])]; 
               x_Rx_mat = [x_Rx_mat;...
                           reshape(x_Tx_wnoise_wdm(idx,:),num_block_interleave,[])];
            end
            Lblc_size_temp = size(L_Rx_mat,2)/nC;
            xblc_size_temp = size(x_Rx_mat,2)/nC;
            for idx = 1:nC
                L_Rx_wdm(idx,:) = reshape(L_Rx_mat(:,(idx-1)*Lblc_size_temp+1:idx*Lblc_size_temp),1,[]);
                x_Rx_wdm(idx,:) = reshape(x_Rx_mat(:,(idx-1)*xblc_size_temp+1:idx*xblc_size_temp),1,[]);
            end
        end 
        nEr1 = [];
        nEr2 = [];
        % channels
        for jdx = 1:nC
            r = L_Rx_wdm(jdx,:); 
            msgcap = [];       
            for idx = 1:Nblock/nC
                if Enable_List_decoding
                    msgcap_block = SSCList_Polar_Decoder(r((idx-1)*N+1:idx*N),N,Q1,F,NL,crcg);
                else
                    msgcap_block = SSC_Polar_Decoder(r((idx-1)*N+1:idx*N),N,Q1,F);
                end
                msgcap = [msgcap msgcap_block];
            end
            C_rx = rx_PAM4_Decode(x_Rx_wdm(jdx,:));
            C_rx2 = ~C_rx;
            nEr1 = [nEr1, sum(abs(Cin_mat(jdx,:)-(C_rx2)))]; 
            nEr2 = [nEr2, sum(abs(msg_mat(jdx,:)-msgcap))];
        end
        
        nEr1 = [nEr1, mean(nEr1)];
        nEr2 = [nEr2, mean(nEr2)];
        % display results
        disp('BER w/o polar code:');
        nEr1/size(Cin_mat,2)
        disp('BER with polar code:');
        nEr2/size(msg_mat,2)
end
