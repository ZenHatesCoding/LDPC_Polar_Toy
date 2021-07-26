% (N,k) polar code encoder with cyclic redundancy checking
% https://freevideolectures.com/course/4202/nptel-ldpc-polar-codes-in-g-standard/28
% written by Zhenping Xing
% zhenping.xing@mail.mcgill.ca

function u = Polar_Encoder_wCRC(msg,N,A,Q1,crcg)
%     crcL = 11;
%     crcg = fliplr([1 1 1 0 0 0 1 0 0 0 0 1]); % CRC polynomial
    crcL = length(crcg)-1;
    K = A + crcL;
    [quot, rem] = gfdeconv([zeros(1,crcL) fliplr(msg)], crcg);
    
    msgcrc = [msg fliplr([rem zeros(1,crcL-length(rem))])];
    
    n = log2(N);
    u = zeros(1,N);
    u(Q1(N-K+1:end))=msgcrc;

    m = 1; % number of bits combined
    for d = n-1:-1:0
        for idx = 1:2*m:N
            a = u(idx:idx+m-1);
            b = u(idx+m:idx+2*m-1);
            u(idx:idx+2*m-1) = [mod(a+b,2) b]; 
        end
        m = m * 2;
    end
end
