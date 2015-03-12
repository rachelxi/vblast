% VBLAST Detection

clear all;          % clear workspace
close all;          % close figue windows

Nt = 4;             % send antenna number
Nr = 4;             % receive antenna number

L = 100;            % frame length
SimTimes = 1e3;%10000;   % static count per SNR

EbN0indB = 0:2:30;       % define Eb/N0 range
%EbN0indB = 50;
ModType = 4;             % modulation mode: 1=BPSK, 4=QPSK, 16=16QAM, 64=64QAM
% C/N = SNR = (Eb/No)*(fb/B) => SNR = 10log10(Eb/N0) + 10log10(fb/B)
SNRindB = EbN0indB + 10*log10(log2(ModType));

thisMode   = modem.qammod(ModType);
thisDemod  = modem.qamdemod(thisMode);
% Average Symbol Energy 1 BPSK, 2 QPSK, 5 16QAM, etc...
Es = (mean(thisMode.Constellation .* conj(thisDemod.Constellation))) ;
SNR = zeros(1,length(SNRindB));

dec_zf          = zeros(L,1);
dec_mmse        = zeros(L,1);
dec_zf_sorted   = zeros(L,1);
dec_mmse_sorted = zeros(L,1);
dec_qr_zf       = zeros(L,1);
dec_zf_sqrd     = zeros(L,1);
dec_mmse_sqrd   = zeros(L,1);
   
EB_zf           = zeros(1,length(SNRindB));
EB_mmse         = zeros(1,length(SNRindB));
EB_zf_sorted    = zeros(1,length(SNRindB));
EB_mmse_sorted  = zeros(1,length(SNRindB));
EB_qr_zf        = zeros(1,length(SNRindB));
EB_qr_mmse      = zeros(1,length(SNRindB));
EB_zf_sqrd      = zeros(1,length(SNRindB));
EB_mmse_sqrd    = zeros(1,length(SNRindB));

fid = fopen('vblast_v3.log','a+');
tic
for index = 1:length(SNRindB)
    SNR(index) = 10^(SNRindB(index)/10);
    sigma = sqrt(Es/SNR(index)); 
    % sigma = sqrt(1/SNR(index));
    
    for simcnt= 1:SimTimes
        txMsgBits  = randi([0,1],[1,log2(ModType)*L]);
        txMsgBitsTmp = reshape(txMsgBits,log2(ModType),[]);
        txMsgBitsInt = bi2de(txMsgBitsTmp.','left-msb')';
        txMapped = modulate(thisMode, txMsgBitsInt);
        
        x = reshape(txMapped,Nt,L/Nt);                                    % reshape data to Nt antennas
        % x = sqrt(1/Es)*x;                                                 % normalization 
        AWGN_noise = sqrt(1/2)*sigma*(randn(Nr,L/Nt)+1j*randn(Nr,L/Nt));  % the Nois Matrix
        
        % ======== do detection procedure ========
        r = zeros(Nr,L/Nr);
        for col_idx = 1:L/Nt
            H = sqrt(1/Nt)*sqrt(1/2)*(randn(Nr,Nt) + 1j*randn(Nr,Nt));    % fast fading Rayleigh channel
            r(:,col_idx) = H*x(:,col_idx) + AWGN_noise(:,col_idx);        % received signal
            rsic = r(:,col_idx);
            
             dec_zf((col_idx-1)*Nt+1:col_idx*Nt)   = vblast_zf(rsic,H,ModType);
%            dec_mmse((col_idx-1)*Nt+1:col_idx*Nt) = vblast_mmse(rsic,H,ModType,sigma);
            
            
            dec_qr_zf((col_idx-1)*Nt+1:col_idx*Nt) = qr_zf_sic(rsic,H,ModType);
%             dec_qr_mmse((col_idx-1)*Nt+1:col_idx*Nt) = qr_mmse_sic_v2(rsic,H,ModType,sigma);
            
            dec_zf_sorted((col_idx-1)*Nt+1:col_idx*Nt) = vblast_zf_sorted(rsic,H,ModType);
%            dec_mmse_sorted((col_idx-1)*Nt+1:col_idx*Nt) = vblast_mmse_sorted(rsic,H,ModType,sigma);
     
             dec_zf_sqrd((col_idx-1)*Nt+1:col_idx*Nt) = qr_zf_sic_sorted(rsic,H,ModType);
%            dec_mmse_sqrd((col_idx-1)*Nt+1:col_idx*Nt) = qr_mmse_sic_sorted_v3(rsic,H,ModType,sigma);
            
%             test_dec = qr_mmse_sic_sorted_v2(rsic,H,ModType,sigma);
        end % end of dection loop
         dec_zf_bin = reshape(de2bi(dec_zf,2,'left-msb')',1,[]);
%        det_mmse_bin = reshape(de2bi(dec_mmse,2,'left-msb')',1,[]);
         dec_zf_sorted_bin = reshape(de2bi(dec_zf_sorted,2,'left-msb')',1,[]);
%        dec_mmse_sorted_bin = reshape(de2bi(dec_mmse_sorted,2,'left-msb')',1,[]);
         dec_qr_zf_bin = reshape(de2bi(dec_qr_zf,2,'left-msb')',1,[]);
%        dec_qr_mmse_bin = reshape(de2bi(dec_qr_mmse,2,'left-msb')',1,[]);
         dec_zf_sqrd_bin = reshape(de2bi(dec_zf_sqrd,2,'left-msb')',1,[]);
%        dec_mmse_sqrd_bin = reshape(de2bi(dec_mmse_sqrd,2,'left-msb')',1,[]);
        
         EB_zf(index)      = EB_zf(index) + sum(abs(dec_zf_bin~=txMsgBits));            
%        EB_mmse(index)    = EB_mmse(index) + sum(abs(det_mmse_bin~=txMsgBits));        
         EB_zf_sorted(index)   = EB_zf_sorted(index) + sum(abs(dec_zf_sorted_bin~=txMsgBits));      
%        EB_mmse_sorted(index) = EB_mmse_sorted(index) +sum(abs(dec_mmse_sorted_bin~=txMsgBits));   
         EB_qr_zf(index)   = EB_qr_zf(index) + sum(abs(dec_qr_zf_bin~=txMsgBits));
%        EB_qr_mmse(index) = EB_qr_mmse(index) + sum(abs(dec_qr_mmse_bin~=txMsgBits));
         EB_zf_sqrd(index) = EB_zf_sqrd(index) + sum(abs(dec_zf_sqrd_bin~=txMsgBits));
%        EB_mmse_sqrd(index) = EB_mmse_sqrd(index) + sum(abs(dec_mmse_sqrd_bin~=txMsgBits));
        
    end %end of simcnt loop
    EbN0indB(index)
%     fprintf(fid,'%s\n',datestr(datetime('now')));
%     fprintf(fid,'SNR=%d,EB_zf=%d,EB_mmse=%d,EB_zf_sorted=%d,EB_mmse_sorted=%d,EB_qr_zf=%d,EB_qr_mmse=%d,EB_zf_sqrd=%d,EB_mmse_sqrd=%d\n', ...
%         EbN0indB(index),EB_zf(index),EB_mmse(index),EB_zf_sorted(index),EB_mmse_sorted(index),...
%         EB_qr_zf(index),EB_qr_mmse(index),EB_zf_sqrd(index),EB_mmse_sqrd(index));
end % end of SNR loop
toc
fclose(fid);

TotalBits = ((L*log2(ModType))*SimTimes);
 BER_zf   = EB_zf./TotalBits;
%BER_mmse = EB_mmse./TotalBits;
 BER_zf_sorted  = EB_zf_sorted./TotalBits;
%BER_mmses_sorted = EB_mmse_sorted./TotalBits;
 BER_qr_zf = EB_qr_zf./TotalBits;
%BER_qr_mmse = EB_qr_mmse./TotalBits;
 BER_zf_sqrd = EB_zf_sqrd./TotalBits;
%BER_mmse_sqrd = EB_mmse_sqrd./TotalBits;

% show the results. MMSE
figure(108);
semilogy(EbN0indB,BER_zf,'-ro','LineWidth',2);hold on;
semilogy(EbN0indB,BER_qr_zf,'-b*','LineWidth',2);hold on;
semilogy(EbN0indB,BER_zf_sqrd,'-kv','LineWidth',2);hold on;
semilogy(EbN0indB,BER_zf_sorted,'-mpentagram','LineWidth',2);hold on;
xlabel('Eb/N_{0} in dB');ylabel('BER');
legend('ZF','ZF-QRD','ZF-SQRD','ZF-BLAST');
grid on;
print(108,'-djpeg','ZF_compare.jpeg');


% show the results. ZF
% figure(101);
% semilogy(EbN0indB,BER_zf,'-ro','LineWidth',2);hold on;
% semilogy(EbN0indB,BER_qr_zf,'-b*','LineWidth',2);hold on;
% semilogy(EbN0indB,BER_zf_sqrd,'-kv','LineWidth',2);hold on;
% semilogy(EbN0indB,BER_zf_sorted,'-mpentagram','LineWidth',2);hold on;
% xlabel('Eb/N_{0} in dB');ylabel('BER');
% legend('ZF','ZF-QRD','ZF-SQRD','ZF-BLAST');
% grid on;
% print(101,'-djpeg','MMSE.jpeg');
