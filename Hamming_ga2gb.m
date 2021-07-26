function M = Hamming_ga2gb(a,b)
    M = [];
    for i = a:b
        M = [M,Hamming_g(i)];
    end
end