% written by Zhenping Xing
% zhenping.xing@mail.mcgill.ca
% ref: Chase - 1972 - A class of algorithms for decoding block codes wit.
% algorithm 2
function xout = Hamming_SIHO_Chase_Decoder(xin,H,d_over2,numpos)
    % try num_error = 2;
    
    x_HD = xin<0;
    syndrome = mod(H*x_HD.',2);
    
    xout = x_HD;
    if sum(syndrome) ~= 0
        [~, pos] = mink(abs(xin),numpos);
        T0 = zeros(size(xin));
        minWeight = Inf;
        for i = 1:d_over2
            error_pos = nchoosek(pos,i);
            for j = 1:size(error_pos,1)
                posj = error_pos(j,:);
                T = T0;
                T(posj)=1;
                xnew = mod(x_HD+T,2);
    
                s = mod(H*xnew.',2);
                if sum(s) == 0                
                    Error = mod(x_HD+xnew,2);
                    Weight = abs(xin)*Error.';
                    if Weight < minWeight
                        minWeight = Weight;
                        xout = xnew;
                    end
                end
            end
        end
        
    end

end