% LDPC decoder using the full (simplified) update rule, inner loop
% completely vectorized but still relatively slow
function xh = decode_LDPC(L, H, iterations)
    [row_i, col_i] = find(H);  % get row and column indices
    
    % initialize variable to check node messages with channel output
    VtoC = sparse(row_i, col_i, L(col_i));

    % main iterations
    for it = 1:iterations
        % compute check to variable sum(CtoV,1)node messages
        VtoC_sign = spfun(@sign, VtoC);
        VtoC_abs = spfun(@abs, VtoC);
        
        phiVtoC = spfun(@(x)(log(coth(x/2))), VtoC_abs);
        phiVtoC_sum = sum(phiVtoC,2);
        
        [tri,~,values] = find(VtoC_sign);
        totalsign_VtoC = accumarray(tri,values,[],@prod);
       
        CtoV_abs = spfun(@(x)(log(coth(x/2))), sparse(row_i, col_i,phiVtoC_sum(row_i)) - phiVtoC);
        CtoV_sign = sparse(row_i, col_i, totalsign_VtoC(row_i)) .* VtoC_sign;
        CtoV = CtoV_sign .* CtoV_abs;        
       
        
        % compute variable to check node messages, pretty simple
        CtoV_sum = sum(CtoV,1);
        VtoC = sparse(row_i, col_i, L(col_i) + CtoV_sum(col_i)) - CtoV;        

        % stopping criterion, all parity checks are fulfilled
        L_total = CtoV_sum + L;
        if all(L_total > 0)
            break;
        end       
    end
    
    % binary decision    
    xh = L_total < 0;
end