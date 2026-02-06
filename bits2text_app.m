function [text, symb] = bits2text_app(empfBits, channel)
%BITS2TEXT_APP  Minimal: convert bitstream back to ASCII text.
%
% Inputs:
%   empfBits : received bitstream (row or column vector)
%   channel    : struct with fields
%              .len_cInfoBits   – header length in bits (may be 0)
%              .origTxtLenBits  – true text length in bits (= nChars*8)
%
% Outputs:
%   text : recovered string
%   symb : corresponding uint8 values (ASCII)

    % Ensure row vector
    empfBits = empfBits(:).';

    % --- 1) Remove header bits (control info) ---
    if isfield(channel, 'len_cInfoBits') && channel.len_cInfoBits > 0
        if numel(empfBits) <= channel.len_cInfoBits
            text = "";
            symb = [];
            warning('bits2text_app: received bits shorter than len_cInfoBits.');
            return;
        end
        dataBits = empfBits(channel.len_cInfoBits+1:end);
    else
        dataBits = empfBits;
    end

    % --- 2) Truncate to original text length (if known) ---
    if isfield(channel, 'origTxtLenBits') && channel.origTxtLenBits > 0
        L = channel.origTxtLenBits;
        if numel(dataBits) < L
            warning('bits2text_app: received bits shorter than origTxtLenBits; truncating to received length.');
            L = numel(dataBits);
        end
        dataBits = dataBits(1:L);
    end

    % Ensure length is a multiple of 8
    L8 = floor(numel(dataBits)/8)*8;
    if L8 == 0
        text = "";
        symb = [];
        warning('bits2text_app: effective data bit length is zero.');
        return;
    end
    dataBits = dataBits(1:L8);

    % --- 3) 8 bits -> 1 byte -> text ---
    mBits = reshape(dataBits, 8, []).';          % [Nchar x 8]
    symb  = uint8(bi2de(mBits, 'left-msb'));     % 0..255

    % Map NULL (0) to space to avoid UI issues
    symb(symb == 0) = 32;

    text = char(symb).';
end