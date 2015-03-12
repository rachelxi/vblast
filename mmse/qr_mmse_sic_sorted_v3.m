% function : qr_mmse_sic_sorted
% description : 
% D.Wubben J.Rinas Efficient Algorithm for Detecting Layered Space-Time Codes,p5
%

function dec = qr_mmse_sic_sorted_v3(r,H,ModType,sigma)
    thisMode   = modem.qammod(ModType);
    thisDemod  = modem.qamdemod(ModType);
    % Es = (mean(thisMode.Constellation .* conj(thisDemod.Constellation))) ;
    
    [Nr,Nt] = size(H);
    
    Hext = [H;sigma.*eye(Nt)];      % extended format of Channel Matrix H
    R = zeros(Nt,Nt);
    Q = Hext;
    k = 1:Nt;
    Qtemp = zeros(Nr+Nt,Nt);
    
    % SQRD Alogrithm
    for i = 1:Nt                        % loop column
        Qtemp(:,i:Nt) = Q(:,i:Nt);
        Qtemp(:,1:i-1) = Inf;
        
        [~,kmin] = min(sum(abs(Qtemp).^2,1));
        k(:,[i,kmin]) = k(:,[kmin,i]);  % culumn exchange
        R(:,[i,kmin]) = R(:,[kmin,i]);  % column exchange
        Q(1:Nr+i-1,[i,kmin]) = Q(1:Nr+i-1,[kmin,i]);  % column exchange here is different!!! 
        
        R(i,i) = norm(Q(:,i));
        Q(:,i) = Q(:,i)/R(i,i);
        
        for l = (i+1):Nt
            R(i,l) = Q(:,i)'*Q(:,l);
            Q(:,l) = Q(:,l) - R(i,l)*Q(:,i);
        end
    end
    permutation = k;
    
    % Signal Detection
    dec_p = zeros(1,length(r));
    r = [r;zeros(Nt,1)];        % extended format of r
    
    d   = zeros(1,Nt);
    z   = zeros(1,Nt);
    c   = zeros(1,Nt);
    
    y = Q'*r;
    for k = Nt:-1:1
        for i = (k+1):Nt
            d(k) = d(k) + R(k,i)*c(i);
        end
        z(k) = y(k) - d(k); 
        dec_p(k) = demodulate(thisDemod,z(k)/R(k,k));
        c(k) = modulate(thisMode,dec_p(k));
    end
    
    dec(permutation) = dec_p;
      
end
