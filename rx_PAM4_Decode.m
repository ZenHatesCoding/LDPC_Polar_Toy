function rx_bits = rx_PAM4_Decode(rx_symbols)



n = length(rx_symbols);

rx_bits = zeros(1,2*n);

% PAM4 [0  1  3  2]
%      -3 -1  1  3
% boundaries: 3/sqrt(10) = 0.9487 2/sqrt(10) = 0.6325 1/sqrt(10) = 0.3162 
table = -3:2:3;
norm_coef = sqrt(norm(table)^2/length(table));

for i = 1:n
    if real(rx_symbols(i)) >= 2/norm_coef
        rx_bits(2*(i-1)+1) = 1;
        rx_bits(2*(i-1)+2) = 0;
    elseif real(rx_symbols(i)) >= 0 && real(rx_symbols(i)) < 2/norm_coef
        rx_bits(2*(i-1)+1) = 1;
        rx_bits(2*(i-1)+2) = 1;
    elseif real(rx_symbols(i)) >= -2/norm_coef && real(rx_symbols(i)) < 0
        rx_bits(2*(i-1)+1) = 0;
        rx_bits(2*(i-1)+2) = 1;
    else 
        rx_bits(2*(i-1)+1) = 0;
        rx_bits(2*(i-1)+2) = 0;
    end

end
end

