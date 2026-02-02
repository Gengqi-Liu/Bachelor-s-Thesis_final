function [numError, errRate] = bitFehlerRaten_app(empfData, sendData, datenTyp)
%BITFEHLERRATEN_APP
% Compute bit error count and bit error rate (BER) between sent and received data.
%
% Inputs:
%   empfData : reconstructed data at receiver
%              - image: uint8 / double [H x W x C]
%              - text : char vector / string
%   sendData : original data at transmitter (e.g. app.SendeDatei)
%   datenTyp : 'Bild' | 'Image' | 'Text' (case-insensitive)
%
% Outputs:
%   numError : number of bit errors
%   errRate  : bit error rate (numError / compared bits)

    dt = string(datenTyp);
    dt = lower(strtrim(dt));

    switch dt
        % ================= Text =================
        case "text"
            % Tx: text -> bits
            txText  = char(sendData);
            symbTx  = double(txText);
            mBitsTx = de2bi(symbTx, 8);
            vBitsTx = reshape(mBitsTx.', 1, []);

            % Rx: text -> bits
            rxText  = char(empfData);
            symbRx  = double(rxText);
            mBitsRx = de2bi(symbRx, 8);
            vBitsRx = reshape(mBitsRx.', 1, []);

            % Compare up to common length
            len_txBits = length(vBitsTx);
            len_rxBits = length(vBitsRx);
            min_Len    = min(len_txBits, len_rxBits);

            if min_Len == 0
                numError = 0;
                errRate  = 0;
                warning('bitFehlerRaten_app: text bit length is 0, returning 0 error rate.');
                return;
            end

            [numError, errRate] = biterr(vBitsTx(1:min_Len), vBitsRx(1:min_Len));

        otherwise
            error('bitFehlerRaten_app: unknown datenTyp = "%s", expected "Text".', datenTyp);
    end
end