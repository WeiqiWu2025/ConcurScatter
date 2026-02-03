function [txWaveform,dataSpMapped] = survey_ConcurScatter_funcOFDMSymDerived(dataBits,cfgFormat,varargin)
% wlanWaveformGenerator WLAN waveform generation
%   WAVEFORM = wlanWaveformGenerator(DATA,CFGFORMAT) generates a waveform
%   for a given format configuration and information bits. The generated
%   waveform contains a single packet with no idle time. For OFDM based
%   formats, the data scrambler initial state is 93 and the packet is
%   windowed for spectral controls with a windowing transition time of 1e-7
%   seconds.
%
%   WAVEFORM is a complex Ns-by-Nt matrix containing the generated
%   waveform, where Ns is the number of time domain samples, and Nt is the
%   number of transmit antennas.
%
%   DATA is the information bits including any MAC padding to be coded
%   across the number of packets to generate, i.e., representing multiple
%   concatenated PSDUs. It can be a double or int8 typed binary vector.
%   Alternatively, it can be a scalar cell array or a vector cell array
%   with length equal to number of users. Each element of the cell array
%   must be a double or int8 typed, binary vector. When DATA is a vector or
%   scalar cell array, it applies to all users. When DATA is a vector cell
%   array, each element applies to a single user. For each user, the bit
%   vector applied is looped if the number of bits required across all
%   packets of the generation exceeds the length of the vector provided.
%   This allows a short pattern to be entered, e.g. [1;0;0;1]. This pattern
%   will be repeated as the input to the PSDU coding across packets and
%   users. The number of data bits taken from a data stream for the ith
%   user when generating a packet is given by the ith element of the
%   CFGFORMAT.PSDULength property times eight.
%
%   CFGFORMAT is a format configuration object of type <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>,
%   <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a>, <a href="matlab:help('wlanHETBConfig')">wlanHETBConfig</a>, <a href="matlab:help('wlanDMGConfig')">wlanDMGConfig</a>, <a href="matlab:help('wlanS1GConfig')">wlanS1GConfig</a>, <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>,
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>. The format of the generated waveform is
%   determined by the type of CFGFORMAT. The properties of CFGFORMAT are
%   used to parameterize the packets generated including the data rate and
%   PSDU length.
%
%   WAVEFORM = wlanWaveformGenerator(DATA,CFGFORMAT,Name,Value) specifies
%   additional name-value pair arguments described below. When a name-value
%   pair is not specified, its default value is used.
%   
%   'NumPackets'               The number of packets to generate. It must
%                              be a positive integer. The default value is
%                              1.
%   
%   'IdleTime'                 The length in seconds of an idle period
%                              after each generated packet. The valid range
%                              depends on the format to generate. For DMG
%                              it must be 0 or greater than or equal to
%                              1e-6 seconds. For all other formats it must
%                              be 0 or greater than or equal to 2e-6
%                              seconds. The default value is 0 seconds.
% 
%   'ScramblerInitialization'  Scrambler initial state(s), applied for HE,
%                              DMG, S1G, VHT, HT, and non-HT OFDM formats.
%                              It must be a double or int8-typed scalar or
%                              matrix containing integer values. The valid
%                              range depends on the format to generate.
%
%                              For DMG Control PHY the valid range is
%                              between 1 and 15 inclusive. For other DMG
%                              formats, the valid range is between 1 and
%                              127 inclusive.
%
%                              For HE, S1G, VHT, and HT formats, the valid
%                              range is between 1 and 127 inclusive.
%
%                              For Non-HT OFDM the valid range depends on
%                              whether bandwidth signaling is enabled.
%
%                              When bandwidth signaling is not used
%                              (CFGFORMAT.SignalChannelBandwidth is false)
%                              the specified value is the initial state of
%                              the scrambler. The specified value must be
%                              between 1 and 127 inclusive.
%
%                              When bandwidth signaling is used
%                              (CFGFORMAT.SignalChannelBandwidth is true),
%                              the specified value is the initial
%                              pseudorandom scrambler sequence as
%                              described in IEEE 802.11-2016 Table 17-7.
%                              The valid range depends on the value of
%                              CFGFORMAT.BandwidthOperation and
%                              CFGFORMAT.ChannelBandwidth. For more
%                              information, see: <a href="matlab:doc('wlanwaveformgenerator')">wlanWaveformGenerator</a>
%                              documentation.
%
%                              To initialize all packets with the same
%                              state for all users, specify this input as a
%                              scalar.
%
%                              To initialize each packet with a distinct
%                              state, specify this input as a column vector
%                              of length NumUsers. The function uses these
%                              initial states for all users.
%
%                              To initialize each packet with a distinct
%                              state for each user, specify this input as a
%                              matrix of size NumPackets-by-NumUsers. Each
%                              column specifies the initial states for a
%                              single user. Each row specifies the initial
%                              state of the corresponding packet.
%
%                              If the number of packets in the waveform
%                              exceeds the number of rows you provide in
%                              this input, the function generates the
%                              waveform by looping the rows.
%
%                              For all formats except DMG, the default
%                              value is 93, which is the example state
%                              given in IEEE Std 802.11-2016 Section
%                              I.1.5.2. For the DMG format, the value
%                              specified will override the
%                              ScramblerInitialization property of the
%                              configuration object. The mapping of the
%                              initialization bits on scrambler schematic
%                              X1 to X7 is specified in IEEE Std
%                              802.11-2016, Section 17.3.5.5. For more
%                              information, see: <a href="matlab:doc('wlanwaveformgenerator')">wlanWaveformGenerator</a>
%                              documentation.
% 
%   'WindowTransitionTime'     The windowing transition length in seconds,
%                              applied to OFDM based formats. For all
%                              formats except DMG it must be a nonnegative
%                              scalar and no greater than 16e-6 seconds.
%                              Specifying it as 0 turns off windowing. For
%                              all formats except DMG, the default value is
%                              1e-7 seconds. For DMG OFDM format it must be
%                              a nonnegative scalar and no greater than
%                              9.6969e-08 (256/2640e6) seconds. The default
%                              value for DMG format is 6.0606e-09
%                              (16/2640e6) seconds.
%   Examples:
%
%   Example 1:
%       %  Generate a time domain signal txWaveform for an 802.11ac VHT
%       %  transmission with 10 packets and 20 microsecond idle period 
%       %  between packets.
%
%       numPkts = 10;                   % 10 packets in the waveform
%      
%       cfgVHT = wlanVHTConfig();       % Create format configuration
%       % Change properties from defaults
%       cfgVHT.NumTransmitAntennas = 2; % 2 transmit antennas
%       cfgVHT.NumSpaceTimeStreams = 2; % 2 spatial streams
%       cfgVHT.MCS = 1;                 % Modulation: QPSK Rate: 1/2 
%       cfgVHT.APEPLength = 1024;       % A-MPDU length in bytes
%
%       % Create bit vector containing concatenated PSDUs
%       numBits = cfgVHT.PSDULength*8*numPkts;
%       dataBits = randi([0 1],numBits,1);
%
%       txWaveform = wlanWaveformGenerator(dataBits, cfgVHT, ...
%           'NumPackets', numPkts, 'IdleTime', 20e-6, ...
%           'WindowTransitionTime', 1e-7);
%
%   Example 2:
%       %  Produce a waveform containing a single 802.11a packet without 
%       %  any windowing.
%          
%       cfgNonHT = wlanNonHTConfig(); % Create format configuration
%
%       psdu = randi([0 1], cfgNonHT.PSDULength*8, 1); % Create a PSDU
%
%       txWaveform = wlanWaveformGenerator(psdu, cfgNonHT, ...
%           'WindowTransitionTime', 0); % Disable windowing
%
%   Example 3:
%       %  Produce a waveform containing a single DMG packet with a
%       %  specified scrambler initialization.
%          
%       cfgDMG = wlanDMGConfig(); % Create format configuration
%       cfgDMG.MCS = 1;           % Single carrier
%       cfgDMG.ScramblerInitialization = 93; % Specify initialization
%
%       psdu = randi([0 1], cfgDMG.PSDULength*8, 1); % Create a PSDU
%
%       txWaveform = wlanWaveformGenerator(psdu, cfgDMG);
%
%   Example 4:
%       %  Produce a waveform containing a multiple DMG packets, each with
%       %  a random scrambler initialization.
%          
%       cfgDMG = wlanDMGConfig(); % Create format configuration
%       numPkts = 4; % Generate 4 packets
%
%       % Create bit vector containing concatenated PSDUs
%       numBits = cfgDMG.PSDULength*8*numPkts;
%       dataBits = randi([0 1],numBits,1);
%
%       txWaveform = wlanWaveformGenerator(dataBits, cfgDMG, ...
%           'NumPackets', numPkts, ...
%           'IdleTime', 1e-5, ...
%           'ScramblerInitialization', randi([1 15],numPkts,1));
%
%   Example 5:
%       %  Produce a waveform containing an 802.11ax HE single user packet 
%       %  without any windowing.
%          
%       cfgHESU = wlanHESUConfig(); % Create format configuration
%
%       psdu = randi([0 1], getPSDULength(cfgHESU)*8, 1); % Create a PSDU
%
%       txWaveform = wlanWaveformGenerator(psdu, cfgHESU, ...
%           'WindowTransitionTime', 0); % Disable windowing
%
%   Example 6:
%       %  Produce a waveform containing an 802.11ax HE multi user packet
%       %  with Packet Extension, for two RUs and two users.
%          
%       cfgHEMU = wlanHEMUConfig([192 192]); % Create format configuration
%       cfgHEMU.User{1}.NominalPacketPadding = 16;
%       cfgHEMU.User{2}.NominalPacketPadding = 8;
%
%       % Generate a random PSDU for each user
%       psdu = cell(1,numel(cfgHEMU.User));
%       psduLength = getPSDULength(cfgHEMU);
%       for i = 1:numel(cfgHEMU.User)
%           psdu{i} = randi([0 1],psduLength(i)*8,1,'int8');
%       end
% 
%       txWaveform = wlanWaveformGenerator(psdu, cfgHEMU);
%
%   Example 7:
%       %  Produce a waveform containing an 802.11ax HE trigger-based 
%       %  packet without any windowing. The trigger method used to 
%       %  generate the HE TB PPDU is set to TriggerFrame (by default).
%
%       cfgHETB = wlanHETBConfig(); % Create format configuration
%
%       psdu = randi([0 1], getPSDULength(cfgHETB)*8, 1); % Create a PSDU
%
%       txWaveform = wlanWaveformGenerator(psdu, cfgHETB, ...
%           'WindowTransitionTime', 0); % Disable windowing
%
%   See also wlanVHTConfig, wlanHTConfig, wlanNonHTConfig, wlanS1GConfig,
%   wlanDMGConfig, wlanHESUConfig, wlanHEMUConfig, wlanHETBConfig,
%   wirelessWaveformGenerator.

%   Copyright 2015-2020 The MathWorks, Inc.

%#codegen

% Check number of input arguments
coder.internal.errorIf(mod(nargin, 2) == 1, 'wlan:wlanWaveformGenerator:InvalidNumInputs');

% Validate the format configuration object is a valid type
validateattributes(cfgFormat,{'wlanVHTConfig','wlanHTConfig','wlanNonHTConfig','wlanS1GConfig','wlanDMGConfig','wlanHESUConfig','wlanHEMUConfig','wlanHETBConfig'},{'scalar'},mfilename,'format configuration object');
s = validateConfig(cfgFormat);

% Get format
isNonHT = isa(cfgFormat,'wlanNonHTConfig');
isDMG = isa(cfgFormat,'wlanDMGConfig');
inDSSSMode = isNonHT && strcmpi(cfgFormat.Modulation,'DSSS');
isDMGOFDM = isDMG && strcmp(phyType(cfgFormat),'OFDM');
isS1G = isa(cfgFormat,'wlanS1GConfig');
isVHT = isa(cfgFormat,'wlanVHTConfig');
isHT = isa(cfgFormat,'wlanHTConfig');
isHEMU = isa(cfgFormat,'wlanHEMUConfig');
isHE = isa(cfgFormat,'wlanHESUConfig') || isa(cfgFormat,'wlanHETBConfig') || isHEMU;
inOFDMMode = isHE || isDMGOFDM || isS1G || isVHT || isHT || ...
    (isNonHT && strcmpi(cfgFormat.Modulation,'OFDM')); 

overrideObjectScramInit = false;

% P-V pairs
% Define default values for WindowTransitionTime
if isDMG
    winTransitTime = 16/2640e6; % Windowing length of 16
else
    winTransitTime = 1e-7;
end

defaultScramblerInitialization = 93;
if isa(cfgFormat,'wlanNonHTConfig') && cfgFormat.SignalChannelBandwidth
    % If bandwidth signaling is used then take only the most significant
    % bits required from the default
    numScramBits = 7;
    [~,numRandomBits] = scramblerRange(cfgFormat);
    defaultScramblerInitialization = bitshift(93,-(numScramBits-numRandomBits));
end

% Default values
defaultParams = struct('NumPackets', 1, ...
                    'IdleTime', 0, ...
                    'ScramblerInitialization', defaultScramblerInitialization, ...
                    'WindowTransitionTime', winTransitTime);

if nargin==2
    useParams = defaultParams;
else          
    % Extract each P-V pair
    if isempty(coder.target) % Simulation path
        p = inputParser;

        % Get values for the P-V pair or set defaults for the optional arguments
        addParameter(p,'NumPackets',defaultParams.NumPackets);
        addParameter(p,'IdleTime',defaultParams.IdleTime);
        addParameter(p,'ScramblerInitialization',defaultParams.ScramblerInitialization);
        addParameter(p,'WindowTransitionTime',defaultParams.WindowTransitionTime);
        % Parse inputs
        parse(p,varargin{:});

        useParams = p.Results;
    else % Codegen path
        pvPairs = struct('NumPackets', uint32(0), ...
                         'IdleTime', uint32(0), ...
                         'ScramblerInitialization', uint32(0), ...
                         'WindowTransitionTime', uint32(0));

        % Select parsing options
        popts = struct('PartialMatching', true);

        % Parse inputs
        pStruct = coder.internal.parseParameterInputs(pvPairs,popts,varargin{:});

        % Get values for the P-V pair or set defaults for the optional arguments
        useParams = struct;
        useParams.NumPackets = coder.internal.getParameterValue(pStruct.NumPackets,defaultParams.NumPackets,varargin{:});
        useParams.IdleTime = coder.internal.getParameterValue(pStruct.IdleTime,defaultParams.IdleTime,varargin{:});
        useParams.ScramblerInitialization = coder.internal.getParameterValue(pStruct.ScramblerInitialization,defaultParams.ScramblerInitialization,varargin{:});
        useParams.WindowTransitionTime = coder.internal.getParameterValue(pStruct.WindowTransitionTime,defaultParams.WindowTransitionTime,varargin{:});
    end

    % Validate each P-V pair
    % Validate numPackets
    validateattributes(useParams.NumPackets,{'numeric'},{'scalar','integer','>=',0},mfilename,'''NumPackets'' value');
    % Validate idleTime
    validateattributes(useParams.IdleTime,{'numeric'},{'scalar','real','>=',0},mfilename,'''IdleTime'' value');
    if isDMG
        minIdleTime = 1e-6;
    else % S1G, VHT, HT, non-HT
        minIdleTime = 2e-6;
    end
    coder.internal.errorIf((useParams.IdleTime > 0) && (useParams.IdleTime < minIdleTime),'wlan:wlanWaveformGenerator:InvalidIdleTimeValue',sprintf('%1.0d',minIdleTime));
    % Validate scramblerInit
    if ~inDSSSMode 
        if isDMG && any(useParams.ScramblerInitialization~=93)
            if strcmp(phyType(cfgFormat),'Control')
                coder.internal.errorIf(any((useParams.ScramblerInitialization<1) | (useParams.ScramblerInitialization>15)),'wlan:wlanWaveformGenerator:InvalidScramblerInitialization','Control',1,15);
            elseif wlan.internal.isDMGExtendedMCS(cfgFormat.MCS)
                % At least one of the initialization bits must be
                % non-zero, therefore determine if the pseudorandom
                % part can be 0 given the extended MCS and PSDU length.
                if all(wlan.internal.dmgExtendedMCSScramblerBits(cfgFormat)==0)
                    minScramblerInit = 1; % Pseudorandom bits cannot be all zero
                else
                    minScramblerInit = 0; % Pseudorandom bits can be all zero
                end
                coder.internal.errorIf(any((useParams.ScramblerInitialization<minScramblerInit) | (useParams.ScramblerInitialization>31)),'wlan:wlanWaveformGenerator:InvalidScramblerInitialization','SC extended MCS',minScramblerInit,31);
            else
                coder.internal.errorIf(any((useParams.ScramblerInitialization<1) | (useParams.ScramblerInitialization>127)),'wlan:wlanWaveformGenerator:InvalidScramblerInitialization','SC/OFDM',1,127);
            end
            overrideObjectScramInit = true;
        else
            if isNonHT && inOFDMMode && ...
              any(strcmp(cfgFormat.ChannelBandwidth,{'CBW20','CBW40','CBW80','CBW160'})) ...
              && cfgFormat.SignalChannelBandwidth
                % Non-HT may include bandwidth signaling
            
                % Validate type
                validateattributes(useParams.ScramblerInitialization,{'double','int8'}, ...
                  {'real','integer','2d','nonempty'},mfilename,'''ScramblerInitialization'' value');
              
                % Validate range
                range = scramblerRange(cfgFormat);
                minVal = range(1);
                maxVal = range(2);
                % Check for correct range
                if any((useParams.ScramblerInitialization<minVal) | (useParams.ScramblerInitialization>maxVal),'all')
                    coder.internal.error('wlan:wlanWaveformGenerator:InvalidScramInitBWSignaling',minVal,maxVal);
                end
            else
                % Validate scrambler initialization
                validateattributes(useParams.ScramblerInitialization,{'double','int8'},{'real','integer','2d','nonempty'},mfilename,'''ScramblerInitialization'' value');
                if any((useParams.ScramblerInitialization<1) | (useParams.ScramblerInitialization>127),'all')
                    coder.internal.error('wlan:wlanWaveformGenerator:InvalidScramInit',1,127);
                end
            end
        end
    end
    % Validate WindowTransitionTime
    if inOFDMMode 
        % Set maximum limits for windowing transition time based on bandwidth and format
        if isHE
            maxWinTransitTime = 6.4e-6; % Seconds
        elseif isDMG
            maxWinTransitTime = 9.6969e-08; % Seconds
        elseif isS1G
            maxWinTransitTime = 16e-6; % Seconds
        elseif isNonHT
            switch cfgFormat.ChannelBandwidth
                case 'CBW5'
                    maxWinTransitTime = 6.4e-6; % Seconds
                case 'CBW10'
                    maxWinTransitTime = 3.2e-6; % Seconds
                otherwise % 'CBW20'
                    maxWinTransitTime = 1.6e-6; % Seconds
            end
        else % HT/VHT
            maxWinTransitTime = 1.6e-6; % Seconds
        end
        validateattributes(useParams.WindowTransitionTime,{'numeric'},{'real','scalar','>=',0,'<=',maxWinTransitTime},mfilename,'''WindowTransitionTime'' value');
    end 
end
windowing = inOFDMMode && useParams.WindowTransitionTime > 0;

if isVHT || isS1G
    numUsers = cfgFormat.NumUsers;  
elseif isHEMU
    allocInfo = ruInfo(cfgFormat);
    numUsers = allocInfo.NumUsers;
else
    numUsers = 1;
end

% Cross validation
coder.internal.errorIf(all(size(useParams.ScramblerInitialization,2) ~= [1 numUsers]),'wlan:wlanWaveformGenerator:ScramInitNotMatchNumUsers');

psduLength = s.PSDULength;

% Validate that data bits are present if PSDULength is nonzero
if iscell(dataBits) % SU and MU
    % Data must be a scalar cell or a vector cell of length Nu
    coder.internal.errorIf(~isvector(dataBits) || all(length(dataBits) ~= [1 numUsers]), 'wlan:wlanWaveformGenerator:InvalidDataCell');
    
    for u = 1:length(dataBits)
        if ~isempty(dataBits{u}) && (psduLength(u)>0) % Data packet
            validateattributes(dataBits{u},{'double','int8'},{'real','integer','vector','binary'},mfilename,'each element in cell data input');
        else
            % Empty data check if not NDP
            coder.internal.errorIf((psduLength(u)>0) && isempty(dataBits{u}),'wlan:wlanWaveformGenerator:NoData');
        end
    end
    if isscalar(dataBits) 
        % Columnize and expand to a [1 Nu] cell
        dataCell = repmat({int8(dataBits{1}(:))},1,numUsers);
    else % Columnize each element
        numUsers = numel(dataBits); % One cell element per user
        dataCell = repmat({int8(1)},1,numUsers); 
        for u = 1:numUsers                
            dataCell{u} = int8(dataBits{u}(:));
        end
    end
else % SU and MU: Data must be a vector
    if ~isempty(dataBits) && any(psduLength > 0) % Data packet
        validateattributes(dataBits,{'double','int8'},{'real','integer','vector','binary'}, mfilename,'Data input');

        % Columnize and expand to a [1 Nu] cell
        dataCell = repmat({int8(dataBits(:))}, 1, numUsers);
    else % NDP
        % Empty data check if not NDP
        coder.internal.errorIf(any(psduLength > 0) && isempty(dataBits),'wlan:wlanWaveformGenerator:NoData');

        dataCell = {int8(dataBits(:))};
    end
end

% Number of bits in a PSDU for a single packet (convert bytes to bits)
numPSDUBits = psduLength*8;

% Repeat to provide initial state(s) for all users and packets
scramInit = repmat(useParams.ScramblerInitialization,1,numUsers/size(useParams.ScramblerInitialization,2)); % For all users
pktScramInit = scramInit(mod((0:useParams.NumPackets-1).',size(scramInit,1))+1, :);

% Get the sampling rate of the waveform
if isDMG
    if strcmp(phyType(cfgFormat),'OFDM')
        sr = 2640e6;
    else
        sr = 1760e6;
    end
    numTxAnt = 1;
    numPktSamples = s.NumPPDUSamples;
elseif isHE
    [psps,trc] = wlan.internal.hePacketSamplesPerSymbol(cfgFormat);
    numPktSamples = psps.NumPacketSamples;
    numTxAnt = cfgFormat.NumTransmitAntennas;
    cbw = wlan.internal.cbwStr2Num(cfgFormat.ChannelBandwidth);
    sr = cbw*1e6;
    giType = cfgFormat.GuardInterval;
    sf = cbw*1e-3; % Scaling factor to convert bandwidth and time in ns to samples
    FFTLen = trc.TDFTHE*sf;
    symLength = trc.TSYM*sf;
    cpLen = trc.TGIData*sf;
    heltfSymLen = trc.THELTFSYM*sf; 
    Npe = wlan.internal.heNumPacketExtensionSamples(trc.TPE,cbw);
elseif inDSSSMode % DSSS format
    sr = 11e6;
    numTxAnt = 1;
    giType = ''; % For codegen
    FFTLen = 0;  % For codegen
    info = wlan.internal.dsssInfo(cfgFormat);
    numPktSamples = info.NumPPDUSamples;

    lstf = []; % For codegen
    lltf = []; % For codegen
    lsig = []; % For codegen
else % NonHT/VHT/HT/S1G OFDM format
    chanBW = cfgFormat.ChannelBandwidth;
    sr = wlan.internal.cbwStr2Num(chanBW)*1e6;
    if isNonHT
        giType = 'Long'; % Always
        FFTLen = wlan.internal.cbw2nfft(chanBW);
    elseif isS1G
        giType = cfgFormat.GuardInterval;
        FFTLen = (sr/2e6)*64;
    else % For VHT/HT formats
        giType = cfgFormat.GuardInterval;
        FFTLen = wlan.internal.cbw2nfft(chanBW);
    end
    if any(strcmp(chanBW,{'CBW10','CBW5'}))
        numTxAnt = 1; % Override and set to 1 only, for 802.11j/p
    else
        numTxAnt = cfgFormat.NumTransmitAntennas;
    end
    numPktSamples = real(s.NumPPDUSamples); % real for codegen

    % Generate the legacy preamble fields for applicable formats
    if ~isS1G
        lstf = wlanLSTF(cfgFormat);
        lltf = wlanLLTF(cfgFormat);
        lsig = wlanLSIG(cfgFormat);
    end
end

if isVHT
    vhtsiga = wlanVHTSIGA(cfgFormat);
    vhtstf = wlanVHTSTF(cfgFormat);
    vhtltf = wlanVHTLTF(cfgFormat);
    vhtsigb = wlanVHTSIGB(cfgFormat);
    preamble = [lstf; lltf; lsig; vhtsiga; vhtstf; vhtltf; vhtsigb];
elseif isHT
    htSig = wlanHTSIG(cfgFormat);
    htstf = wlanHTSTF(cfgFormat);
    htltf = wlanHTLTF(cfgFormat);
    preamble = [lstf; lltf; lsig; htSig; htstf; htltf];
elseif isNonHT
    if strcmp(cfgFormat.Modulation,'OFDM')
        preamble = [lstf; lltf; lsig];
    else % DSSS
        preamble = [wlan.internal.wlanDSSSPreamble(cfgFormat); wlan.internal.wlanDSSSHeader(cfgFormat)];
    end
elseif isS1G
    if ~strcmp(packetFormat(cfgFormat),'S1G-Long')
        stf = wlan.internal.s1gSTF(cfgFormat);
        [ltf1, ltf2n] = wlan.internal.s1gLTF(cfgFormat);
        sig = wlan.internal.s1gSIG(cfgFormat);
        preamble = [stf; ltf1; sig; ltf2n];
    else % Preamble == 'Long'
        stf = wlan.internal.s1gSTF(cfgFormat);
        ltf1 = wlan.internal.s1gLTF1(cfgFormat);
        siga = wlan.internal.s1gSIGA(cfgFormat);
        dstf = wlan.internal.s1gDSTF(cfgFormat);
        dltf = wlan.internal.s1gDLTF(cfgFormat);
        sigb = wlan.internal.s1gSIGB(cfgFormat);
        preamble = [stf; ltf1; siga; dstf; dltf; sigb];
    end
elseif isDMG
    if strcmp(phyType(cfgFormat),'OFDM')
        % In OFDM PHY preamble fields are resampled to OFDM rate
        preamble = wlan.internal.dmgResample([wlan.internal.dmgSTF(cfgFormat); wlan.internal.dmgCE(cfgFormat)]);
        brfields = wlan.internal.dmgResample([wlan.internal.dmgAGC(cfgFormat); wlan.internal.dmgTRN(cfgFormat)]);
    else
        preamble = [wlan.internal.dmgSTF(cfgFormat); wlan.internal.dmgCE(cfgFormat)];
        brfields = [wlan.internal.dmgAGC(cfgFormat); wlan.internal.dmgTRN(cfgFormat)];
    end
elseif isHE
    LSTF = wlan.internal.heLSTF(cfgFormat);
    LLTF = wlan.internal.heLLTF(cfgFormat);
    LSIG = wlan.internal.heLSIG(cfgFormat);
    RLSIG = LSIG; % RL-SIG is identical to L-SIG
    SIGA = wlan.internal.heSIGA(cfgFormat);
    STF = wlan.internal.heSTF(cfgFormat);
    LTF = wlan.internal.heLTF(cfgFormat);
    if isHEMU
        SIGB = wlan.internal.heSIGB(cfgFormat);
        preamble = [LSTF; LLTF; LSIG; RLSIG; SIGA; SIGB; STF; LTF];
    else
        preamble = [LSTF; LLTF; LSIG; RLSIG; SIGA; STF; LTF];
    end
end

if windowing
    % Calculate parameters for windowing
    wlength = 2*ceil(useParams.WindowTransitionTime*sr/2);
    bLen = wlength/2; % Number of samples overlap at the end of the packet
    if isDMG
        % No waveform extension for the non-OFDM fields in the preamble
        aLen = 0;
        % No waveform extension due to windowing when BRP field is present
        bLen = bLen*(~wlan.internal.isBRPPacket(cfgFormat));
        windowedPktLength = numPktSamples+bLen;
    else
        aLen = bLen-1; % Number of samples overlap at start of packet
        windowedPktLength = numPktSamples+wlength-1;
    end
else    
    % Define unused windowing variables for codegen
    wlength = 0;
    windowedPktLength = numPktSamples+wlength-1;
    aLen = 0;
    bLen = 0;
end

% Define a matrix of total simulation length
numIdleSamples = round(sr*useParams.IdleTime);
pktWithIdleLength = numPktSamples+numIdleSamples;
txWaveform = complex(zeros(useParams.NumPackets*pktWithIdleLength,numTxAnt));

for i = 1:useParams.NumPackets
    % Extract PSDU for the current packet
    psdu = getPSDUForCurrentPacket(dataCell, numPSDUBits, i);
    
    % Generate the PSDU with the correct scrambler initial state
    if isVHT
        if any(cfgFormat.APEPLength > 0)
            data = wlanVHTData(psdu,cfgFormat,pktScramInit(i,:));
        else % NDP
            data = complex(zeros(0,cfgFormat.NumTransmitAntennas));
        end
    elseif isHT
        if cfgFormat.PSDULength > 0                    
            [data,dataSpMapped] = funcWlanHTData(psdu{1},cfgFormat,pktScramInit(i,:));
        else % NDP or sounding packet
            data = complex(zeros(0,cfgFormat.NumTransmitAntennas));
        end
    elseif isNonHT
        if strcmp(cfgFormat.Modulation, 'OFDM')
            data = wlanNonHTData(psdu{1},cfgFormat,pktScramInit(i,:));
        else % DSSS
            data = wlan.internal.wlanDSSSData(psdu{1},cfgFormat);
        end
    elseif isS1G
        data = wlan.internal.s1gData(psdu,cfgFormat,pktScramInit(i,:));
    elseif isDMG
        % Header and data scrambled so generate for each packet together
        
        % Override scrambler initialization in configuration object if
        % supplied by the user to the waveform generator
        if overrideObjectScramInit
            cfgFormat.ScramblerInitialization = pktScramInit(i,:);
        end
        
        data = [wlan.internal.dmgHeader(psdu{1},cfgFormat); wlan.internal.dmgData(psdu{1},cfgFormat); brfields];
    elseif isHE
        if any(psduLength > 0)
            data = wlan.internal.heData(psdu,cfgFormat,pktScramInit(i,:));
            
            % Midamble processing
            Mma = cfgFormat.MidamblePeriodicity;
            Nsym = s.NumDataSymbols;
            Nma = wlan.internal.numMidamblePeriods(cfgFormat,Nsym); % Midamble period
            if Nma>0
                % Reshape data symbols till last midamble in to data symbol blocks 
                dataSymBlk = reshape(data(1:Nma*symLength*Mma,:),symLength*Mma,Nma,numTxAnt);
                % Repeat HELTF symbols for each data symbol block
                heltfSymBlk = repmat(permute(LTF,[1,3,2]),1,Nma,1);
                % Append midamble after each data symbol block
                dataMidambleBlk = [dataSymBlk; heltfSymBlk];
                % Reshape and append leftover data samples after the last midamble
                dataMidambleBlkReshape = permute(reshape(dataMidambleBlk,[],1,numTxAnt),[1 3 2]);
                data = [dataMidambleBlkReshape(:,:,1); data(Nma*symLength*Mma+1:end,:)]; % Index 3rd dimension for codegen
            end
            dataPacket = data;
            
            % Packet Extension
            lastDataSymBlk = data(end-symLength+cpLen+1:end,:);
            packetExt = getPacketExtensionData(lastDataSymBlk,Npe);
            data = [dataPacket; packetExt];
        else % NDP
            lastDataSymBlk = preamble(end-heltfSymLen+cpLen+1:end,:);
            packetExt = getPacketExtensionData(lastDataSymBlk,Npe);
            data = packetExt;
        end
    end

    % Construct packet from preamble and data
    packet = [preamble; data];

    if windowing
        % Window each packet
        if isDMG
            windowedPacket = wlan.internal.dmgWindowing(packet,wlength,cfgFormat);
        elseif isHE
            windowedPacket = windowingFunction(packet,psps.NumSamplesPerSymbol,psps.CPPerSymbol,wlength,numTxAnt);
        else
            if isS1G && strcmp(cfgFormat.GuardInterval,'Short')
                % For S1G the first data symbol is always Long GI
                numSamplesSymbol = 40e-6*sr; %40 us from Table 24.4 IEEE P802.11ah/D5.0
                numSamplesBeforeShortGI = size(preamble,1)+numSamplesSymbol;
            else
                % For other formats short GI begins after the preamble
                numSamplesBeforeShortGI = size(preamble,1);
            end
            windowedPacket = wlan.internal.wlanWindowing(packet,FFTLen,wlength,giType,numSamplesBeforeShortGI);
        end

        % Overlap-add the windowed packets
        if useParams.NumPackets==1 && numIdleSamples==0 % Only one packet which wraps     
            txWaveform = windowedPacket(aLen+(1:numPktSamples), :);
            % Overlap start of packet with end
            tmp = windowedPacket(end-bLen+1:end,:);
            txWaveform(1:bLen,:) = txWaveform(1:bLen,:)+tmp(1:bLen,:);
            % Overlap end of packet with start
            txWaveform(end-aLen+1:end,:) = txWaveform(end-aLen+1:end,:)+windowedPacket(1:aLen,:);
        else
            if i==1 % First packet (which wraps)
                % First packet wraps to end of waveform
                txWaveform(1:(numPktSamples+bLen),:) = windowedPacket(1+aLen:end,:);
                txWaveform(end-aLen+1:end,:) = windowedPacket(1:aLen,:);
            elseif i==useParams.NumPackets && numIdleSamples==0 % Last packet which wraps
                % Last packet wraps to start of waveform
                startIdx = (i-1)*pktWithIdleLength-aLen+1;
                txWaveform(startIdx:end,:) = txWaveform(startIdx:end,:)+windowedPacket(1:end-bLen,:);
                tmp = windowedPacket(end-bLen+1:end,:);
                txWaveform(1:bLen,:) = txWaveform(1:bLen,:)+tmp(1:bLen,:);
            else % Packet does not wrap
                % Normal windowing overlap between packets
                idx = (i-1)*pktWithIdleLength-aLen+(1:windowedPktLength);
                txWaveform(idx,:) = txWaveform(idx,:)+windowedPacket;
            end
       end
    else
        % Construct entire waveform
        txWaveform((i-1)*pktWithIdleLength+(1:numPktSamples),:) = packet;
    end
end
end

function packetExt = getPacketExtensionData(lastDataSymBlk,Npe)
    % Cyclic extension of last symbol for packet extension
    if size(lastDataSymBlk,1)>=Npe
        packetExt = lastDataSymBlk(1:Npe,:);
    else
        buffCntt = ceil(Npe/size(lastDataSymBlk,1));
        dataBuffer = repmat(lastDataSymBlk,buffCntt,1);
        packetExt = dataBuffer(1:Npe,:);
    end
end

function psdu = getPSDUForCurrentPacket(dataCell,numPSDUBitsPerPacket,packetIdx)
    numUsers = length(dataCell); % == length(numPSDUBits)
    psdu = repmat({int8(1)},1,numUsers); % Cannot use cell(1, numUsers) for codegen
    for u = 1:numUsers
        psdu{u} = wlan.internal.parseInputBits(dataCell{u},numPSDUBitsPerPacket(u),(packetIdx-1)*numPSDUBitsPerPacket(u));
    end
end

function y = windowingFunction(x,samplesPerSymbol,cpPerSymbol,wLength,Nt)
    % windowingFunction(...) returns the time-domain windowed signal for
    % the OFDM signal. The windowing function for OFDM waveform is defined
    % in IEEE Std 802.11-2016.
    assert(size(x,1)==sum(samplesPerSymbol));
    assert(all(size(cpPerSymbol)==size(samplesPerSymbol)));

    % Window length must be less than or equal to twice the CP length (ignore zeros)
    coder.internal.errorIf(wLength>(2*min(cpPerSymbol(cpPerSymbol>0))), ...
        'wlan:wlanWindowing:InvalidWindowLength');

    Ns = size(x,1); % Number of samples
    Nsym = size(samplesPerSymbol,2); % Number of OFDM symbols

    % Allocate output, the extra samples allow for rampup and down
    y = complex(zeros(Ns+wLength-1,Nt));

    % Offset in samples of each OFDM symbol
    startOffset = cumsum([0 samplesPerSymbol]);

    % For each OFDM symbol extract the portions which overlap, cyclic extend
    % and apply windowing equation. Preallocate additional first and last to
    % create rampup and rampdown
    prefixOverlap = complex(zeros(wLength-1,Nsym+2,Nt));
    postfixOverlap = complex(zeros(wLength-1,Nsym+2,Nt));
    if coder.target('MATLAB')
        for i = 1:Nsym
            % Standard defined windowing equation for each sample in extended
            % symbol
            [~,w] = wlan.internal.windowingEquation(wLength,samplesPerSymbol(i));

            % Extract data symbol
            dataSym = x(startOffset(i)+(1:samplesPerSymbol(i)),:);

            % Extend the symbol with a prefix to create the desired window
            % transition and apply windowing equation
            prefixOverlap(:,i+1,:) = permute( ...
                [dataSym((1:(wLength/2-1))+(samplesPerSymbol(i)-cpPerSymbol(i)-(wLength/2-1)),:); ...
                dataSym(1:(wLength/2),:)], ...
                [1 3 2]).*w(1:(wLength-1));

            % Extend the symbol with a postfix to create the desired window
            % transition and apply windowing equation
            postfixOverlap(:,i+1,:) = permute( ...
                [dataSym(end-(wLength/2-1)+(1:(wLength/2-1)),:); ...
                dataSym(cpPerSymbol(i)+(1:wLength/2),:)], ...
                [1 3 2]).*w(end-(wLength-2):end);
        end
    else
        for i = 1:Nsym
            % Standard defined windowing equation for each sample in extended
            % symbol
            [~,w] = wlan.internal.windowingEquation(wLength,samplesPerSymbol(i));

            % Extract data symbol
            dataSym = x(startOffset(i)+(1:samplesPerSymbol(i)),:);

            % Extend the symbol with a prefix to create the desired window
            % transition and apply windowing equation
            for j = 1:Nt
                prefixOverlap(:,i+1,j) = permute( ...
                    [dataSym((1:(wLength/2-1))+(samplesPerSymbol(i)-cpPerSymbol(i)-(wLength/2-1)),j); ...
                    dataSym(1:(wLength/2),j)], ...
                    [1 3 2]).*w(1:(wLength-1));

                % Extend the symbol with a postfix to create the desired window
                % transition and apply windowing equation
                tmp = w(end-(wLength-2):end);
                postfixOverlap(:,i+1,j) = permute( ...
                    [dataSym(end-(wLength/2-1)+(1:(wLength/2-1)),j); ...
                    dataSym(cpPerSymbol(i)+(1:wLength/2),j)], ...
                    [1 3 2]).*tmp(1:(wLength-2+1));
            end
        end
    end

    % Overlap the prefix and postfix regions, note first prefix region
    overlap = prefixOverlap(:,2:end,:)+postfixOverlap(:,1:end-1,:);

    % First samples at output will be the rampup i.e. overlap with zeros
    y(1:wLength/2-1,:) = overlap(1:wLength/2-1,1,:);

    % Construct windowed symbols from overlap regions and symbol samples
    for i = 1:Nsym
        % Extract symbol from input
        dataSym = x(startOffset(i)+(1:samplesPerSymbol(i)),:);

        % Extract start, middle and end portions and store
        startPortion = permute(overlap(wLength/2:end,i,:),[1 3 2]);
        middlePortion = dataSym((wLength/2)+1:end-(wLength/2-1),:);
        endPortion = permute(overlap(1:(wLength/2-1),i+1,:),[1 3 2]);
        idx = wLength/2-1+startOffset(i)+(1:samplesPerSymbol(i));
        y(idx,:) = [startPortion; middlePortion; endPortion];
    end

    % Last samples output will be rampdown i.e. overlap with zeros
    idx = wLength/2-1+startOffset(Nsym+1)+(1:wLength/2);
    y(idx,:) = permute(overlap(wLength/2:end,Nsym+1,:),[1 3 2]);
end