function [mEmpfDataBits, data] = EigenModeRx_app(mFrameRxNoCP, params, channel, len_cInfoBits)
% EigenModeRx_app
%
% EigenMode receiver (text-only) with:
%   - Channel estimation: Zero Forcing (ZF) only
%   - Equalization: MMSE only
%   - Demapping: BPSK header + QAM payload
%
% Inputs:
%   mFrameRxNoCP  [iNfft x iNewNoBlocks x iNoRxAnt]
%   params        struct, must contain V1 [iNoTxAnt x iNoTxAnt x iNfft]
%   channel         struct (kept for interface compatibility; modes enforced here)
%   len_cInfoBits header length in bits (BPSK)
%
% Outputs:
%   mEmpfDataBits recovered bits (row vector)
%   data          debug struct

%% Enforce fixed modes (project constraint)
channel.estimator = 'Zero Forcing';
channel.equalizer = 'MMSE';

%% 1) Unpack parameters
iNoBlocks    = params.iNoBlocks;
iNfft        = params.iNfft;
iModOrd      = params.iModOrd;
iNoTxAnt     = params.iNoTxAnt;
iNoRxAnt     = params.iNoRxAnt;
iNoSubBlocks = params.iNoSubBlocks;

if ~isfield(params,'V1') || isempty(params.V1)
    error('EigenModeRx_app: V1 (precoder) missing in params.');
end
V1 = params.V1; % [Tx x Tx x SC]

mFrameRx = mFrameRxNoCP;

%% 2) Split preambles and data blocks
AnzSubFrames = floor(iNoBlocks / iNoSubBlocks);

vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
vPreambleFreq = (1/sqrt(iNfft)) * fft(vPreambleTime, iNfft);

% Preamble blocks: [SC x SubFrames x Rx x Tx]
mPreambleRx = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);
for tx = 1:iNoTxAnt
    mPreambleRx(:,:,:,tx) = ...
        mFrameRx(1:iNfft, ...
                 tx : (iNoSubBlocks+iNoTxAnt) : ...
                 (iNoSubBlocks*(AnzSubFrames-1)+AnzSubFrames*iNoTxAnt), :);
end

% Data blocks: [SC x Blocks x Rx]
mDataRx = zeros(iNfft, iNoBlocks, iNoRxAnt);
v = 0;
for sf = 1:AnzSubFrames
    if sf == AnzSubFrames
        blkRange = iNoSubBlocks*(sf-1)+1 : iNoSubBlocks*sf + rem(iNoBlocks,iNoSubBlocks);
        idxRange = iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v+rem(iNoBlocks,iNoSubBlocks);
    else
        blkRange = iNoSubBlocks*(sf-1)+1 : iNoSubBlocks*sf;
        idxRange = iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v;
    end
    mDataRx(:, blkRange, :) = mFrameRx(1:iNfft, idxRange, :);
    v = v + iNoTxAnt + iNoSubBlocks;
end

%% 3) SNR estimate (from trailing blocks)
if size(mFrameRx,2) >= 10
    mZeroRx = mFrameRx(:, end-9:end, :);
else
    mZeroRx = mFrameRx;
end

mPreambleRxFreq = (1/sqrt(iNfft)) * fft(mPreambleRx, iNfft, 1);

% Noise power per (SC,Rx)
mPowNoise = squeeze(mean(abs(mZeroRx).^2, 2));  % [SC x Rx]

% Signal power per (SC,Rx,Tx) from last preamble of each Tx
mPowSig = squeeze(abs(mPreambleRxFreq(:,end,:,:)).^2); % [SC x Rx x Tx]

mSNR = zeros(iNfft, iNoRxAnt, iNoTxAnt);
for tx = 1:iNoTxAnt
    for rx = 1:iNoRxAnt
        mSNR(:,rx,tx) = (mPowSig(:,rx,tx) - mPowNoise(:,rx)) ./ max(mPowNoise(:,rx), eps);
    end
end
mSNR(~isfinite(mSNR) | mSNR<=0) = 1e4;
iSNR = mean(mSNR(:));

%% 4) Channel estimation (ZF only)
mCTF = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);
mCIR = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);

for tx = 1:iNoTxAnt
    for rx = 1:iNoRxAnt
        for sf = 1:AnzSubFrames
            rxCol = mPreambleRxFreq(:,sf,rx,tx);
            H = rxCol ./ vPreambleFreq;
            mCTF(:,sf,rx,tx) = H;
            mCIR(:,sf,rx,tx) = sqrt(iNfft) * ifft(H);
        end
    end
end

%% 5) Equalization (MMSE only) on effective channel Heff = H * V1
mDataRxFreq = (1/sqrt(iNfft)) * fft(mDataRx, iNfft, 1); % [SC x Blocks x Rx]
mDataRxEq   = zeros(iNfft, iNoBlocks, iNoTxAnt);

for sf = 1:AnzSubFrames
    if sf == AnzSubFrames
        blkRange = (sf-1)*iNoSubBlocks+1 : sf*iNoSubBlocks + rem(iNoBlocks,iNoSubBlocks);
    else
        blkRange = (sf-1)*iNoSubBlocks+1 : sf*iNoSubBlocks;
    end

    for b = blkRange
        for sc = 1:iNfft
            H  = squeeze(mCTF(sc,sf,:,:));      % [Rx x Tx]
            He = H * V1(:,:,sc);                % [Rx x Tx]
            y  = squeeze(mDataRxFreq(sc,b,:));  % [Rx x 1]

            A = He' * He + eye(iNoTxAnt) * (1./max(iSNR, eps));
            if any(~isfinite(A(:))) || rcond(A) < 1e-10
                A = He' * He + eye(iNoTxAnt) * 1e-3;
            end

            x = A \ (He' * y);
            mDataRxEq(sc,b,:) = sqrt(iNoTxAnt) * x;
        end
    end
end

%% 6) Symbols to bits (text-only: BPSK header + QAM payload)
vSym = reshape(mDataRxEq, 1, []);
vSym(~isfinite(vSym)) = 0;

if len_cInfoBits > 0
    vHdrBits = qamdemod(vSym(1:len_cInfoBits), 2, 'UnitAveragePower', true);
else
    vHdrBits = [];
end

vPayloadSym = vSym(len_cInfoBits+1:end);
idxPay      = qamdemod(vPayloadSym, 2^iModOrd, 'UnitAveragePower', true);
bitsPayMat  = de2bi(idxPay, iModOrd, 'right-msb').';
vBitsPay    = bitsPayMat(:);

mEmpfDataBits = [vHdrBits(:); vBitsPay(:)].';

%% 7) Debug output
data = struct();
data.channelUeb     = mCTF;
data.channelimpuls  = mCIR;
data.mDataRxEq    = permute(mDataRxEq, [3 1 2]); % [Tx x SC x Block]
data.iSNR         = iSNR;
data.AnzSubFrames = AnzSubFrames;
end