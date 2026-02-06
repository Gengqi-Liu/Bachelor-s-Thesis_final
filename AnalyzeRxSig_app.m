function [mEmpfDataBits, data] = AnalyzeRxSig_app(rxFrame, params, channel)
%ANALYZERXSIG_APP (text-only)
% RX top-level:
%   1) Downconvert to baseband
%   2) Schmidl & Cox timing
%   3) Cut OFDM blocks and remove CP
%   4) Dispatch to MIMO-mode-specific RX

    % ---- Enforce text-only + fixed estimator/equalizer ----
    channel.estimator  = 'Zero Forcing';
    channel.equalizer  = 'MMSE';

    % ---- Unpack parameters ----
    iNfft        = params.iNfft;
    iNg          = params.iNg;
    iNb          = params.iNb;
    iModOrd      = params.iModOrd;
    iNoRxAnt     = params.iNoRxAnt;
    iNewNoBlocks = params.iNewNoBlocks;

    fBBFreq      = params.fBBFreq;
    fDACFreq     = params.fDACFreq;
    fCarrFreq    = params.fCarrFreq;

    mimoModeStr  = string(params.mimoMode);

    if isfield(channel,'len_cInfoBits') && ~isempty(channel.len_cInfoBits)
        len_cInfoBits = channel.len_cInfoBits;
    else
        len_cInfoBits = 0;
    end

    mEmpfDataBits = [];
    data          = struct();

    %% 2) Downconvert to baseband
    mRecFrame = rxFrame; % [Nr x Nsamp]

    vBB1  = DeModulateSignal_app(mRecFrame(1,:), fBBFreq, fDACFreq, fCarrFreq);
    Lbb   = numel(vBB1);

    mDataRxDem = zeros(iNoRxAnt, Lbb);
    mDataRxDem(1,:) = vBB1;

    for iC = 2:iNoRxAnt
        vBB = DeModulateSignal_app(mRecFrame(iC,:), fBBFreq, fDACFreq, fCarrFreq);
        if numel(vBB) ~= Lbb
            L = min(Lbb, numel(vBB));
            mDataRxDem(iC,1:L) = vBB(1:L);
        else
            mDataRxDem(iC,:) = vBB;
        end
    end

    %% 3) Schmidl & Cox timing
    vFrStart   = zeros(1, iNoRxAnt);
    vCarOff    = zeros(1, iNoRxAnt);
    MetricData = cell(1, iNoRxAnt);

    for iC = 1:iNoRxAnt
        [frameStartTmp, fCarOffTmp, metricStruct] = ...
            EstFrameStart_app(mDataRxDem(iC,:), iNfft, iNg);
        vFrStart(iC)   = frameStartTmp;
        vCarOff(iC)    = fCarOffTmp;
        MetricData{iC} = metricStruct;
    end

    if all(~isfinite(vFrStart)) || isempty(vFrStart)
        iFrSync = 10;
    else
        iFrSync = min(vFrStart(isfinite(vFrStart)));
    end

    safetyMargin = 10;
    iFrStart = iFrSync + iNfft - safetyMargin;
    iFrStart = max(1, min(iFrStart, Lbb));

    % n = 1:Lbb;
    % for iC = 1:iNoRxAnt
    %     cfo = vCarOff(iC);
    %     if isfinite(cfo) && ~isnan(cfo) && abs(cfo) > 0
    %         mDataRxDem(iC,:) = mDataRxDem(iC,:) .* exp(-1j*2*pi/iNfft * cfo .* n);
    %     end
    % end

    %% 4) Cut OFDM blocks and remove CP
    maxSamplesAvail = Lbb - iFrStart + 1;
    maxBlocksAvail  = floor(maxSamplesAvail / iNb);

    if maxBlocksAvail <= 0
        error('AnalyzeRxSig_app: not enough samples for one OFDM block.');
    end

    if maxBlocksAvail < iNewNoBlocks
        iNewNoBlocks = maxBlocksAvail;
    end

    totalNeeded = iNewNoBlocks * iNb;

    mFrameRxTp = zeros(totalNeeded, iNoRxAnt);
    mFrameRx   = zeros(iNb, iNewNoBlocks, iNoRxAnt);

    for iC = 1:iNoRxAnt
        seg = mDataRxDem(iC, iFrStart : iFrStart + totalNeeded - 1);
        mFrameRxTp(:,iC) = seg.';
        mFrameRx(:,:,iC) = reshape(seg, iNb, []);
    end

    mFrameRxNoCP = mFrameRx;
    mFrameRxNoCP(1:iNg,:,:) = [];

    %% 5) Dispatch to mode-specific RX
    switch lower(mimoModeStr)
        case 'alamouti'
            [mEmpfDataBits, dataMode] = AlamoutiRx_app(mFrameRxNoCP, params, channel, len_cInfoBits);

        case {'spatial multiplexing','v-blast','vblast'}
            [mEmpfDataBits, dataMode] = SM_VBLAST_Rx_app(mFrameRxNoCP, params, channel, len_cInfoBits);

        case 'eigenmode'
            [mEmpfDataBits, dataMode] = EigenModeRx_app(mFrameRxNoCP, params, channel, len_cInfoBits);

        otherwise
            error('AnalyzeRxSig_app: unsupported MIMO mode: %s', mimoModeStr);
    end

    %% 6) Fill common debug info
    data = dataMode;
    data.mRecFrame    = mRecFrame;
    data.fDACFreq     = fDACFreq;
    data.fBBFreq      = fBBFreq;
    data.fCarrFreq    = fCarrFreq;
    data.iNfft        = iNfft;
    data.iModOrd      = iModOrd;
    data.ta_TP        = 1 / fBBFreq;
    data.MetricData   = MetricData;
    data.iNewNoBlocks = iNewNoBlocks;
    data.iNg          = iNg;
    data.iNoTxAnt     = params.iNoTxAnt;
    data.iNoRxAnt     = params.iNoRxAnt;
    data.iFrStartUsed = iFrStart;
end


%% ========================================================================
%% Subfunction: SM / V-BLAST RX (preamble extraction + channel est + MMSE EQ)
%% ========================================================================
function [mEmpfDataBits, data] = SM_VBLAST_Rx_app(mFrameRxNoCP, params, channel, len_cInfoBits)
% mFrameRxNoCP : [iNfft x iNewNoBlocks x Nr]

    % ---- Enforce fixed modes ----
    channel.estimator = 'Zero Forcing';
    channel.equalizer = 'MMSE';

    iNoBlocks     = params.iNoBlocks;
    iNfft         = params.iNfft;
    iModOrd       = params.iModOrd;
    iNoTxAnt      = params.iNoTxAnt;
    iNoRxAnt      = params.iNoRxAnt;
    iNoSubBlocks  = params.iNoSubBlocks;

    mimoModeStr   = string(params.mimoMode);

    AnzSubFrames = floor(iNoBlocks / iNoSubBlocks);

    vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
    vPreambleFreq = 1/sqrt(iNfft) * fft(vPreambleTime, iNfft);

    % ---- Extract preambles ----
    mPreambleRx = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);
    for iCTxAnt = 1:iNoTxAnt
        mPreambleRx(:,:,:,iCTxAnt) = ...
            mFrameRxNoCP(1:iNfft, ...
                         iCTxAnt : iNoSubBlocks+iNoTxAnt : iNoSubBlocks*(AnzSubFrames-1)+AnzSubFrames*iNoTxAnt, ...
                         :);
    end

    % ---- Extract data blocks ----
    v = 0;
    mDataRx = [];
    for ii = 1:AnzSubFrames
        if ii == AnzSubFrames
            mDataRx(:, iNoSubBlocks*(ii-1)+1 : iNoSubBlocks*ii + rem(iNoBlocks,iNoSubBlocks), :) = ...
                mFrameRxNoCP(1:iNfft, ...
                             iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v+rem(iNoBlocks,iNoSubBlocks), ...
                             :);
        else
            mDataRx(:, iNoSubBlocks*(ii-1)+1 : iNoSubBlocks*ii, :) = ...
                mFrameRxNoCP(1:iNfft, ...
                             iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v, ...
                             :);
        end
        v = v + iNoTxAnt + iNoSubBlocks;
    end

    % ---- Noise blocks for SNR estimation ----
    if size(mFrameRxNoCP,2) >= 10
        mZeroRx = mFrameRxNoCP(:, end-10+1:end, :);
    else
        mZeroRx = mFrameRxNoCP;
    end

    mPreambleRxFreq = 1/sqrt(iNfft) * fft(mPreambleRx, iNfft, 1);

    mPowNoise    = zeros(iNfft, iNoRxAnt);
    mPowSigNoisy = zeros(iNfft, iNoRxAnt, iNoTxAnt);
    mSNR         = zeros(iNfft, iNoRxAnt, iNoTxAnt);

    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            for k = 1:iNfft
                mPowNoise(k,iCRxAnt) = mean(abs(mZeroRx(k,:,iCRxAnt)).^2);
            end
            mPowSigNoisy(:,iCRxAnt,iCTxAnt) = abs(mPreambleRxFreq(:,end,iCRxAnt,iCTxAnt)).^2;
            mSNR(:,iCRxAnt,iCTxAnt) = ...
                (mPowSigNoisy(:,iCRxAnt,iCTxAnt) - mPowNoise(:,iCRxAnt)) ./ ...
                 max(mPowNoise(:,iCRxAnt), eps);
        end
    end

    badMask = ~isfinite(mSNR) | (mSNR <= 0);
    mSNR(badMask) = 1;

    vSNRperSCperBl = zeros(iNfft,1);
    for k = 1:iNfft
        vSNRperSCperBl(k) = mean(mean(squeeze(mSNR(k,:,:))));
    end

    iSNR = mean(vSNRperSCperBl(:));
    if ~isfinite(iSNR) || iSNR <= 0
        iSNR = 1e4;
    end

    % ---- Channel estimation: ZF only ----
    vPreambleFreqCol = vPreambleFreq(:);
    mCTF = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);
    mCIR = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);

    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            for iCB = 1:AnzSubFrames
                rxCol = mPreambleRxFreq(:,iCB,iCRxAnt,iCTxAnt);
                mCTF(:,iCB,iCRxAnt,iCTxAnt) = rxCol ./ vPreambleFreqCol;
                mCIR(:,iCB,iCRxAnt,iCTxAnt) = sqrt(iNfft) * ifft(mCTF(:,iCB,iCRxAnt,iCTxAnt));
            end
        end
    end

    % ---- Equalization ----
    mDataRxEq = zeros(iNoTxAnt, iNfft, iNoBlocks);

    if strcmpi(mimoModeStr,'v-blast') || strcmpi(mimoModeStr,'vblast')
        channel.equalizer = 'MMSE';
        mDataRxEq = VBlastRx_app( ...
            mDataRx, ...
            mCTF, ...
            struct( ...
                'iNoBlocks',    iNoBlocks, ...
                'iNfft',        iNfft, ...
                'iNoTxAnt',     iNoTxAnt, ...
                'iNoRxAnt',     iNoRxAnt, ...
                'iNoSubBlocks', iNoSubBlocks, ...
                'iModOrd',      iModOrd ...
            ), ...
            iSNR, ...
            channel);
    else
        mDataRxFreq = 1/sqrt(iNfft) * fft(mDataRx, iNfft, 1);
        mDataRxFreq = permute(mDataRxFreq, [3,1,2]); % [Rx x SC x Block]

        for iSF = 1:AnzSubFrames
            if iSF == AnzSubFrames
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks + rem(iNoBlocks,iNoSubBlocks);
            else
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks;
            end

            for iCB = iCB_range
                for k = 1:iNfft
                    mH = squeeze(mCTF(k,iSF,:,:)); % [Rx x Tx]
                    A  = mH' * mH + eye(iNoTxAnt) * (1./iSNR);

                    if any(~isfinite(A(:))) || rcond(A) < 1e-10
                        A = mH' * mH + eye(iNoTxAnt) * 1e-3;
                    end

                    mDataRxEq(:,k,iCB) = sqrt(iNoTxAnt) * ( A \ (mH' * mDataRxFreq(:,k,iCB)) );
                end
            end
        end
    end

    % ---- Optional CPE correction (decision-directed) ----
    M = 2^iModOrd;
    alpha = 0.65;
    thetaHat = zeros(iNoTxAnt,1);

    for iCB = 1:iNoBlocks
        for tx = 1:iNoTxAnt
            s = squeeze(mDataRxEq(tx,:,iCB)).';
            s(~isfinite(s)) = 0;

            mag = abs(s);
            good = mag > 1e-3;
            if nnz(good) < 16
                continue;
            end

            sGood = s(good);
            magGood = mag(good);

            thr = median(magGood);
            keep = magGood >= thr;
            sUse = sGood(keep);

            if numel(sUse) < 16
                continue;
            end

            idx  = qamdemod(sUse, M, 'UnitAveragePower', true);
            sRef = qammod(idx, M, 'UnitAveragePower', true);

            ph = angle(sum(sUse .* conj(sRef)));
            thetaHat(tx) = alpha*thetaHat(tx) + (1-alpha)*ph;

            mDataRxEq(tx,:,iCB) = mDataRxEq(tx,:,iCB) .* exp(-1j*thetaHat(tx));
        end
    end

    % ---- Symbols -> bits (text-only: BPSK header + QAM payload) ----
    mDataRxDet = zeros(iNfft, iNoBlocks, iNoTxAnt);
    for iCTxAnt = 1:iNoTxAnt
        mDataRxDet(:,:,iCTxAnt) = squeeze(mDataRxEq(iCTxAnt,:,:));
    end

    vDataRxDet = reshape(mDataRxDet, 1, iNfft*iNoBlocks*iNoTxAnt);
    vDataRxDet(~isfinite(vDataRxDet)) = 0;

    if len_cInfoBits > 0
        vHeaderBits = qamdemod(vDataRxDet(1:len_cInfoBits), 2, 'UnitAveragePower', true);
    else
        vHeaderBits = [];
    end

    if len_cInfoBits < numel(vDataRxDet)
        vPayloadSym = qamdemod(vDataRxDet(len_cInfoBits+1:end), 2^iModOrd, 'UnitAveragePower', true);
        mBitsPay    = de2bi(vPayloadSym, iModOrd, 'right-msb').';
        vBitsPay    = mBitsPay(:);
    else
        vBitsPay = [];
    end

    mEmpfDataBits = [vHeaderBits(:); vBitsPay(:)].';

    % ---- Debug ----
    data = struct();
    data.channelUeb     = mCTF;
    data.channelimpuls  = mCIR;
    data.mDataRxEq    = mDataRxEq;
    data.AnzSubFrames = AnzSubFrames;
end