function [bits, AnzUeb, metaTxt, original_data_bits] = text2bits_app(SendeDatei, codingMode, anzMaxBits, iModOrd)
    %TEXT2BITS_APP
    % Text to bitstream with optional convolutional coding.
    
    codingMode = string(codingMode);
    
    % 1) Text -> bits
    symb       = double(char(string(SendeDatei)));
    mBits      = de2bi(symb, 8, 'left-msb');
    len_txBits = numel(mBits);
    tx_vBits   = reshape(mBits.', 1, len_txBits);
    
    original_data_bits = tx_vBits;
    metaTxt               = struct();
    metaTxt.txt_size      = len_txBits;
    
    % 2) Control info (payload length)
    tx_infoBits = de2bi(len_txBits, 64, 'left-msb');
    
    trel_ci  = poly2trellis(3, [6 7]);
    tblen_ci = 4;
    
    tx_infoBits_with_tail = [tx_infoBits, randi([0 1], 1, tblen_ci)];
    tx_infoBits_code      = convenc(tx_infoBits_with_tail, trel_ci);
    
    metaTxt.len_cInfoBits = length(tx_infoBits_code);
    
    % 3) Payload coding
    vBits_to_encode = tx_vBits;
    
    switch codingMode
        case "None"
            code_vBits = tx_vBits;
        otherwise
            warning('Unknown coding mode: %s. Using uncoded.', codingMode);
            code_vBits = tx_vBits;
    end
    
    % 4) Prepend control info (repeat iModOrd times)
    for i = 1:iModOrd
        code_vBits = [tx_infoBits_code, code_vBits]; %#ok<AGROW>
    end
    metaTxt.lenCodeBits = length(code_vBits);
    
    % 5) Frame count and padding
    if ~(isscalar(anzMaxBits) && isnumeric(anzMaxBits) && anzMaxBits > 0)
        error('anzMaxBits must be a positive scalar.');
    end
    
    AnzUeb           = ceil(metaTxt.lenCodeBits / anzMaxBits);
    len_total_needed = AnzUeb * anzMaxBits;
    padding_bits     = len_total_needed - metaTxt.lenCodeBits;
    
    if padding_bits > 0
        bits = [code_vBits, randi([0 1], 1, padding_bits)];
    else
        bits = code_vBits;
    end
end