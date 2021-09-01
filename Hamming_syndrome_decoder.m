function xout = Hamming_syndrome_decoder(xin,H)
    s = mod(H*(xin<0).',2);
    xout = xin<0;
    [~,I] = ismember(H.',s.','rows');
    kdx = find(I==1);
    if kdx 
        xout(kdx) = ~xout(kdx);
    end
end