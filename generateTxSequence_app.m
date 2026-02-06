function [mFrameTxCar, meta] = generateTxSequence_app(params, source)
%GENERATETXSEQUENCE_APP
% Main Tx signal generation (text-only):
%   1) text -> bits -> coding -> QAM
%   2) overwrite leading symbols with BPSK header (legacy design)
%   3) Chu preamble
%   4) per-MIMO-mode Tx mapping

    rng(0,'twister');

    % ---- Unpack parameters ----
    iNoBlocks   = params.iNoBlocks;
    iNfft       = params.iNfft;
    iNg         = params.iNg; %#ok<NASGU>
    iNb         = params.iNb; %#ok<NASGU>
    iModOrd     = params.iModOrd;
    iNoTxAnt    = params.iNoTxAnt;

    mimoModeStr = string(params.mimoMode);
    coding      = string(params.codingMode);

    datenTyp = lower(string(source.DatenTyp));
    if datenTyp ~= "text"
        error('Text-only build: source.DatenTyp must be "Text".');
    end

    metaBits = struct();

    % ---- Max bits per frame ----
    anzMaxBits = iModOrd * iNoTxAnt * iNoBlocks * iNfft;

    % ---- Text -> bits (with channel coding) ----
    [vInfoBits, AnzUeb, metaTxt, original_data_bits] = ...
        text2bits_app(string(source.SendeDatei), coding, anzMaxBits, iModOrd);

    len_cInfoBits = metaTxt.len_cInfoBits;

    metaBits.AnzUeb           = AnzUeb;
    metaBits.len_cInfoBits    = len_cInfoBits;
    metaBits.origTxtLenBits   = metaTxt.txt_size;
    metaBits.lenCodeBits      = metaTxt.lenCodeBits;
    metaBits.originalDataBits = original_data_bits;
    metaBits.anzMaxBits       = anzMaxBits;
    metaBits.DatenTyp         = 'text';

    % ---- Use exactly one frame worth of bits ----
    if numel(vInfoBits) < anzMaxBits
        vBits = zeros(1, anzMaxBits);
        vBits(1:numel(vInfoBits)) = vInfoBits;
    else
        vBits = vInfoBits(1:anzMaxBits);
    end

    % ---- Serial-to-parallel and symbol indexing ----
    % vBits -> mBits: [iModOrd x (iNoBlocks*iNfft*iNoTxAnt)]
    mBits = reshape(vBits, iModOrd, iNoBlocks*iNfft*iNoTxAnt);

    % Keep legacy text behavior
    msbOrder = 'right-msb';

    mSymb = bi2de(mBits.', msbOrder); % column vector indices

    % Frequency-domain symbol indices: [iNfft x iNoBlocks x iNoTxAnt]
    mDataTxFreqTp = reshape(mSymb, iNfft, iNoBlocks, iNoTxAnt);

    % ---- QAM modulation ----
    M = 2^iModOrd;
    mDataTxFreq = zeros(size(mDataTxFreqTp));
    for iCAnt = 1:iNoTxAnt
        mDataTxFreq(:,:,iCAnt) = qammod(mDataTxFreqTp(:,:,iCAnt), M, ...
                                        'UnitAveragePower', true);
    end

    % ---- Control information (header) overwrite with BPSK ----
    if len_cInfoBits > 0
        vHeaderBits   = vInfoBits(1:len_cInfoBits);
        vHeaderSymb   = qammod(vHeaderBits, 2, 'UnitAveragePower', true); % BPSK

        lenFrameSymb  = iNfft * iNoBlocks * iNoTxAnt;      % symbols per frame
        vFrameFlat    = reshape(mDataTxFreq, 1, lenFrameSymb);

        if len_cInfoBits <= lenFrameSymb
            vFrameFlat(1:len_cInfoBits) = vHeaderSymb;
        else
            vFrameFlat(:) = vHeaderSymb(1:lenFrameSymb);
        end

        mDataTxFreq = reshape(vFrameFlat, iNfft, iNoBlocks, iNoTxAnt);
    end

    % ---- Chu preamble (time domain) ----
    if exist('ChuSeq','file')
        vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
    else
        vPreambleTime = sqrt(iNoTxAnt) * (randn(iNfft,1) + 1j*randn(iNfft,1))/sqrt(2);
    end

    % ---- Mode-specific Tx mapping ----
    modeLower = lower(strtrim(mimoModeStr));
    switch modeLower
        case 'eigenmode'
            [mFrameTxCar, metaMode] = generateTx_EigenMode_app(params, mDataTxFreq, vPreambleTime);
        case 'alamouti'
            [mFrameTxCar, metaMode] = generateTx_Alamouti_app(params, mDataTxFreq, vPreambleTime);
        case {'spatial multiplexing','v-blast'}
            [mFrameTxCar, metaMode] = generateTx_SM_app(params, mDataTxFreq, vPreambleTime);
        otherwise
            [mFrameTxCar, metaMode] = generateTx_SM_app(params, mDataTxFreq, vPreambleTime);
    end

    meta = mergeMeta(metaBits, metaMode);
end

function metaOut = mergeMeta(metaBase, metaMode)
    metaOut = metaBase;
    if isempty(metaMode)
        return;
    end
    fns = fieldnames(metaMode);
    for k = 1:numel(fns)
        fn = fns{k};
        metaOut.(fn) = metaMode.(fn);
    end
end