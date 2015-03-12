% function vblast_zf.m
% description : lineral zero forcing 
%

function  dec = vblast_zf(r,H,ModType)
   
    thisDemod  = modem.qamdemod(ModType);
    G = (H'*H)\H';      % G = pinv(H)； Moore-Penrose pesudoinverse
    y = G*r;
    dec = demodulate(thisDemod,y);
    dec = dec';
end
% Matrix G : Nt*Nr
