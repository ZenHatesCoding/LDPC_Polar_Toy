% 0 <= i <= 127
function s = Hamming_g(i)
   s = zeros(9,1);
   s(9) = 1;
   for idx = 1:7
       s(idx) = mod(i,2);
       i = floor(i/2);
   end
   sbar = ~s(1:3);
   s(8) = (s(1)&s(3))|(sbar(1)&sbar(2)& sbar(3))|(s(1)& s(2) & sbar(3));
end