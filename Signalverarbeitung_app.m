function data = Signalverarbeitung_app(procParam)
% Signalverarbeitung_app (text-only)
%
% End-to-end RX processing pipeline:
% baseband conversion, sync, channel estimation, equalization, bit recovery,
% text reconstruction, and BER.
%
% Input procParam fields (expected):
%   rxSignal          [Nsamp x Nr]
%   fs               sample rate
%   iNfft, iNg, iNb
%   iNoBlocks, iNewNoBlocks, iNoSubBlocks
%   iNoTxAnt, iNoRxAnt
%   iModOrd
%   fBBFreq, fCarrFreq
%   mimoMode
%   codingMode
%   channelEstimator
%   equalizerMode
%   len_cInfoBits     (optional)
%   origTxtLenBits    (optional)
%   DatenTyp          must be 'Text'
%   SendeDatei        original text payload (optional, for BER)
%
% Output data struct:
%   text
%   numBitError, bitErrorRate
%   plus intermediate results returned by AnalyzeRxSig_app

    % ---- Enforce text-only operation ----
    if ~isfield(procParam,'DatenTyp') || lower(string(procParam.DatenTyp)) ~= "text"
        error('Signalverarbeitung_app: Text-only build. procParam.DatenTyp must be "Text".');
    end

    % ---- Unpack inputs ----
    rxSignal = procParam.rxSignal;   % [Nsamp x Nr]
    fs       = procParam.fs;

    % ---- Build parameter structs for PHY processing ----
    params = struct( ...
        'iNoBlocks',    procParam.iNoBlocks, ...
        'iNfft',        procParam.iNfft, ...
        'iNg',          procParam.iNg, ...
        'iNb',          procParam.iNb, ...
        'iModOrd',      procParam.iModOrd, ...
        'iNoTxAnt',     procParam.iNoTxAnt, ...
        'iNoRxAnt',     procParam.iNoRxAnt, ...
        'iNewNoBlocks', procParam.iNewNoBlocks, ...
        'iNoSubBlocks', procParam.iNoSubBlocks, ...
        'fBBFreq',      procParam.fBBFreq, ...
        'fDACFreq',     fs, ...
        'fCarrFreq',    procParam.fCarrFreq, ...
        'mimoMode',     procParam.mimoMode ...
    );

    if isfield(procParam,'len_cInfoBits') && ~isempty(procParam.len_cInfoBits)
        len_cInfoBits = procParam.len_cInfoBits;
    else
        len_cInfoBits = 0;
    end

    if isfield(procParam,'origTxtLenBits') && ~isempty(procParam.origTxtLenBits)
        origTxtLenBits = procParam.origTxtLenBits;
    else
        origTxtLenBits = [];
    end

    kanal = struct( ...
        'len_cInfoBits',  len_cInfoBits, ...
        'code',           procParam.codingMode, ...
        'schaetzer',      procParam.channelEstimator, ...
        'entzerrer',      procParam.equalizerMode, ...
        'origTxtLenBits', origTxtLenBits ...
    );

    % ---- Pass-through EigenMode matrices if provided ----
    if isfield(params,'mimoMode') && strcmpi(strtrim(lower(string(params.mimoMode))), 'eigenmode')
        if isstruct(procParam) && isfield(procParam,'V1') && ~isempty(procParam.V1)
            params.V1 = procParam.V1;
        end
        if isstruct(procParam) && isfield(procParam,'U1') && ~isempty(procParam.U1)
            params.U1 = procParam.U1;
        end
        if isstruct(procParam) && isfield(procParam,'S1') && ~isempty(procParam.S1)
            params.S1 = procParam.S1;
        end
        if isstruct(procParam) && isfield(procParam,'EigenMeta') && ~isempty(procParam.EigenMeta)
            params.EigenMeta = procParam.EigenMeta;
        end
    end

    % ---- PHY processing: RF -> BB -> sync -> channel -> equalization -> bits ----
    % AnalyzeRxSig_app expects rxFrame as [Nr x Nsamp]
    [mEmpfDataBits, dataChan] = AnalyzeRxSig_app(rxSignal.', params, kanal);

    data = dataChan;

    % ---- Bits -> text reconstruction ----
    rxBits = mEmpfDataBits(:).';     % force row
    L8 = floor(numel(rxBits)/8) * 8; % drop incomplete last byte
    rxBits = rxBits(1:L8);

    [textDecoded, ~] = bits2text_app(rxBits, kanal);
    data.text = textDecoded;

    % ---- BER ----
    if isfield(procParam,'SendeDatei') && ~isempty(procParam.SendeDatei)
        [data.numBitError, data.bitErrorRate] = ...
            bitFehlerRaten_app(textDecoded, procParam.SendeDatei, 'Text');
    else
        data.numBitError  = NaN;
        data.bitErrorRate = NaN;
    end
end