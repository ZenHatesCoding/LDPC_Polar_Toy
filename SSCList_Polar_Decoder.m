% (N,k) polar code  
% simplified successive cancellation list decoder
% ref:
% https://freevideolectures.com/course/4202/nptel-ldpc-polar-codes-in-g-standard/28
% written by Zhenping Xing
% zhenping.xing@mail.mcgill.ca


function msgcap = SSCList_Polar_Decoder(r,N,Q1,F,NL,crcg)
    n = log2(N);
        % depth \in [0,n]
        % node \in [0,2^depth-1]
    
    % storage init
    L = zeros(NL,n+1,N); % beliefs in NL decoders
    Ucap = zeros(NL,n+1,N); % decisions in NL decoders
    PML = Inf*ones(NL,1); % Path Metrics
    PML(1) = 0;
    ns = zeros(1, 2*N-1); % node state
                          % 0 - not activated
                          % 1 - finished L step
                          % 2 - finished R step
                          % 3 - finished U step
    
    f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b)); % minsum
    g = @(a,b,c) b+(1-2*c).*a; % rep code
    
    node = 0; depth = 0;                      
    L(:,1,:) = repmat(r,NL,1,1); % belief of root
    done = 0;
    
    while (done == 0)
        if depth == n
            % leaf
            DM = squeeze(L(:,n+1,node+1)); %decision metrics
            if any(F == (node +1)) % check if frozen
                Ucap(:,n+1,node + 1) = 0; 
                PML = PML + abs(DM).*(DM<0);
            else
                dec = DM < 0;
                PM2 = [PML; PML+abs(DM)];
                % in PM2, first NL are as per DM, next NL are opposite of
                % DM
                
                [PML, pos] = mink(PM2, NL);
                pos1 = pos> NL;
                pos(pos1) = pos(pos1)- NL;
                dec = dec(pos);
                dec(pos1) = 1- dec(pos1);
                L = L(pos,:,:);
                Ucap = Ucap(pos,:,:);
                Ucap(:,n+1,node+1)= dec;
            end
            if node == (N-1) % check if at the last node of the binary tree
                done = 1;
            else
                node = floor(node/2); depth = depth - 1;
            end
        else
            % non leaf
            npos = (2^depth-1)+ node +1; % position of node in node state vec
            switch ns(npos)
                case 0 % Step L: go to left child
                    temp = 2^(n-depth);
                    Ln = squeeze(L(:,depth+1,temp*node+1:temp*(node+1))); % incoming belief
                    a = Ln(:,1:temp/2); b = Ln(:,temp/2+1:end);
                    node = node*2; depth = depth + 1; % next node: left child
                    temp = temp/2;
                    L(:,depth+1, temp*node+1:temp*(node+1)) = f(a,b); % minsum and storage
                    ns(npos) = 1;
                case 1 % Step R: go to right child
                    temp = 2^(n-depth);
                    Ln = squeeze(L(:,depth+1,temp*node+1:temp*(node+1))); % incoming belief
                    a = Ln(:,1:temp/2); b = Ln(:,temp/2+1:end);
                    
                    lnode = node*2; ldepth = depth + 1; % left child
                    ltemp = temp/2;
                    Ucapn = squeeze(Ucap(:,ldepth+1,ltemp*lnode+1:ltemp*(lnode+1))); % incoming decisions from left child
                    
                    node = node*2 + 1; depth = depth + 1; % next node: right child
                    temp = temp/2;
                    L(:,depth+1, temp*node+1:temp*(node+1)) = g(a,b,Ucapn); % minsum and storage
                    ns(npos) = 2;
                case 2 % Step U: go to parent child
                    temp = 2^(n-depth);
                    lnode = node*2; rnode = node*2+1; cdepth = depth + 1; % left and right child    
                    ctemp = temp/2;
                    Ucapl = squeeze(Ucap(:,cdepth+1,ctemp*lnode+1:ctemp*(lnode+1))); % incoming decisions from left child
                    Ucapr = squeeze(Ucap(:,cdepth+1,ctemp*rnode+1:ctemp*(rnode+1))); % incoming decisions from right child
                    Ucap(:,depth+1, temp*node+1:temp*(node+1)) = [mod(Ucapl+Ucapr,2),Ucapr]; % combine
                    node = floor(node/2); depth = depth-1;
                    ns(npos) = 3;
            end
        end
    end
    
    % check CRC
    
    % length(F) = N-K
    msgcapList = squeeze(Ucap(:,n+1,Q1(length(F)+1:end)));
    
    cout = 1; % candidate codeword to be outputted, initially choose the one with the best PM
%     crcg = fliplr([1 1 1 0 0 0 1 0 0 0 0 1]);
    for c1 = 1:NL
        [quot,rem] = gfdeconv(fliplr(msgcapList(c1,:)),crcg);
        if isequal(rem, 0)
            cout = c1;
            break;
        end
    end
    K = N-length(F);
    A = K-length(crcg)+1;
    msgcap = msgcapList(cout,1:A);
end