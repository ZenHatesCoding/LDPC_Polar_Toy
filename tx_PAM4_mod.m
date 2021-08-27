function  mod_symbols = tx_PAM4_mod(bits_in)
   full_len = length(bits_in);
   table = -3:2:3;
   norm_coef = sqrt(norm(table)^2/length(table));
   table = table/norm_coef; 
   table = table([0 1 3 2]+1); % Gray Code mapping pattern

   inp=reshape(bits_in,2,full_len/2); % PAM4 2 bits/symbol
                                      % 2 row full_len/2 col, #col =
                                      % #symbol
   
   mod_symbols=table([2 1]*inp+1);  % maps transmitted bits into PAM4 
end 