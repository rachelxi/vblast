% function qr_zf_sic
% description : zero forcing with QR decomposition.
%

function dec = qr_zf_sic(r,H,ModType)
    thisMode   = modem.qammod(ModType);
    thisDemod  = modem.qamdemod(ModType);
    % Es = (mean(thisMode.Constellation .* conj(thisDemod.Constellation))) ;
    
    [~,Nt] = size(H);
        
    dec = zeros(1,length(r));
    d   = zeros(1,Nt);
    z   = zeros(1,Nt);
    c   = zeros(1,Nt);
    
    [Q,R] = qr(H);
    y = Q'*r;
    
    for k = Nt:-1:1
        for i = (k+1):Nt
            d(k) = d(k) + R(k,i)*c(i);
        end
        z(k) = y(k) - d(k); 
        dec(k) = demodulate(thisDemod,z(k)/R(k,k));
        c(k) = modulate(thisMode,dec(k));
    end
end
