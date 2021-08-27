function dec_symbols = rx_PAM4_Decision(rx_symbols)



n = length(rx_symbols);

dec_symbols = zeros(1,n); 
table = -3:2:3;
norm_coef = sqrt(norm(table)^2/length(table));

for i = 1:n
    if real(rx_symbols(i)) >= 2/norm_coef
        dec_symbols(i) = dec_symbols(i)+3/norm_coef;
    elseif real(rx_symbols(i)) >= 0 && real(rx_symbols(i)) < 2/norm_coef
        dec_symbols(i) = dec_symbols(i)+1/norm_coef;
    elseif real(rx_symbols(i)) >= -2/norm_coef && real(rx_symbols(i)) < 0
        dec_symbols(i) = dec_symbols(i)-1/norm_coef;
    else 
        dec_symbols(i) = dec_symbols(i)-3/norm_coef; 
    end    

end

end

