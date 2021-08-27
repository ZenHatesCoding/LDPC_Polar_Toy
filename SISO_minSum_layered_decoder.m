% reference:
% https://freevideolectures.com/course/4202/nptel-ldpc-polar-codes-in-g-standard/21
% Input:
% r: 1-by-N received coded sequence, r in [-1,1] 
% H, (N-K)-by-N Parity-Checking matrix for linear block code (N,K)
% nlayer: number of layers 
% niter: number of iterations
% offset: offsetting min1 and min2
% Output
% C: 1-by-N coded sequence with error corrected C in [1,0]
% written by Zhenping Xing zhenping.xing@mail.mcgill.ca
function C = SISO_minSum_layered_decoder(r,H,nlayer,niter,Offset)
    if size(r,1)>1
        r = r.';
    end
    C = r;
    L = zeros(size(H));
    num_row_layer = size(H,1)/nlayer; % it has to be an integer

    for i_iter = 1:niter
        for i_layer = 1:nlayer
            idx_row_layer = ((i_layer-1)*num_row_layer+1):...
                            (i_layer*num_row_layer);
            % row operation
            [C,L_layer] = SISO_minSum_row_operation(C,H(idx_row_layer,:),L(idx_row_layer,:),Offset);
            L(idx_row_layer,:) = L_layer;
            % column operation
            [C] = SISO_minSum_col_operation(C,L(idx_row_layer,:));
            
        end

    end    
%     C = (C<0); % mapping from [-1,1] to [1 0];
end

function [C,L] = SISO_minSum_row_operation(C,H,L_old,Offset)
    C = C - sum(L_old,1);
    L = repmat(C,size(H,1),1).*H;
    
    for i = 1:size(L,1)
        % row by row
        Row = L(i,:);
        Row_sign = sign(Row);
        Row_Parity = prod(Row_sign(Row_sign~=0));
        abs_Row = abs(Row);

        sorted_abs_Row = sort(abs_Row(abs_Row>0));
        min1 = sorted_abs_Row(1); 
        min2 = sorted_abs_Row(2);
        % offset
        min1_offset = min1-Offset;
        if min1_offset < 0
            min1_offset = 0;
        end
        
        min2_offset = min2-Offset;
        if min2_offset < 0
            min2_offset = 0;
        end
        Row = min1_offset*sign(Row)*Row_Parity;
        pos_min1 = find(abs_Row == min1);
        Row(pos_min1) = min2_offset*Row_sign(pos_min1)*Row_Parity;
        L(i,:) = Row;
    end
end

function [C] = SISO_minSum_col_operation(C,L_old)    
    C = C + sum(L_old,1);
end