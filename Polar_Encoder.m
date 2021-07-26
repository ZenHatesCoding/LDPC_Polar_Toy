% (N,k) polar code encoder
% https://freevideolectures.com/course/4202/nptel-ldpc-polar-codes-in-g-standard/28
% written by Zhenping Xing
% zhenping.xing@mail.mcgill.ca

function u = Polar_Encoder(msg,N,K,Q1)
    n = log2(N);
    u = zeros(1,N);
    u(Q1(N-K+1:end))=msg;

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
