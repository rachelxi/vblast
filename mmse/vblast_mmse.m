% function vblast_mmse.m
% description :
%

function  dec = vblast_mmse(r,H,ModType,sigma)
    [~,Nt] = size(H);
    thisDemod  = modem.qamdemod(ModType);
    G = (H'*H+sigma.^2*eye(Nt))\H';
    y = G*r;
    dec = demodulate(thisDemod,y);
    dec = dec';
end
% Matrix G : Nt*Nr

