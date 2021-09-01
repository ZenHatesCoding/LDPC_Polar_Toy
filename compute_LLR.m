function L = compute_LLR(y,SE,ro)
    Tx_Table = zeros(1,2^SE);
    switch SE
        case 2 
            for i = 1:2^SE
                Tx_Table(i) = tx_PAM4_mod(fliplr(de2bi(i-1,SE)));
            end
        case 3
            for i = 1:2^SE
                Tx_Table(i) = tx_PAM8_mod(fliplr(de2bi(i-1,SE)));
            end
    end
    N = length(y);
    L = zeros(1,2*N);
    for i = 1:N
        yi = y(i);
        for k = 1:SE
            Lnum = 0;
            Lden = 0;
            for j = 1:2^SE
                bit_j = fliplr(de2bi(j-1,SE));
                if bit_j(k) == 0
                    Lnum = Lnum + exp(-ro/2*(abs(yi-Tx_Table(j))^2));
                else
                    Lden = Lden + exp(-ro/2*(abs(yi-Tx_Table(j))^2));
                end
            end
            L((i-1)*SE+k) = log(Lnum/Lden);
        end
    end
end