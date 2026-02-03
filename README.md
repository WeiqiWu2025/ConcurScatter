# A Review on Concurrent Backscatter Communication over WiFi: Performance Comparison and Future Perspectives

This repository contains the reproduction code for the Paper "Exploiting subcarrier redundancy for concurrent OFDM backscatter communication".

Developed by Weiqi Wu et al. (weiqiwu@ustc.edu.cn), Copyright © 2025 Weiqi Wu et al.

# Code Structure

The code is implemented in a modular manner: the WiFi transmitter first generates excitation signals, which propagate through the transmitter-to-tag channel to the tag. The tag modulates the incident excitation signal, and the modulated signal then propagates through the tag-to-receiver channel to the receiver. The receiver decodes each tag’s data according to the corresponding system design. Each module can be modified to implement customized functionality or to perform further evaluations of the schemes.

Taking file survey_ConcurScatter_Impact_of_backscatter_channel as an example, this file consists of parameter settings, excitation generation, transmitter-to-tag link, tag data generation, tag modulation, tag-to-receiver link, signal superposition from multiple tags, and demodulation.
- **Parameter settings**:
  + Create a format configuration object for a 1-by-1 HT transmission
    ```
    cfgHT = wlanHTConfig;
    cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
    cfgHT.NumTransmitAntennas = 1; % 1 transmit antennas
    cfgHT.NumSpaceTimeStreams = 1; % 1 space-time streams
    cfgHT.PSDULength = 500; % PSDU length in bytes % 64
    cfgHT.MCS = 0; % 1 spatial streams, BPSK rate-1/2
    cfgHT.ChannelCoding = 'BCC'; % BCC channel coding
  + Set SNR, the number of concurrent tags, and the number of packets for each SNR level
    ```
    snr = [0,5,15,25,35];
    global numTags;
    numTags = 3;
    maxNumPackets = 5000;
- **Excitation generation**: Based on parameter settings, the WiFi transmitter generates the corresponding excitation signal.
  ```
  txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
  tx = wlanWaveformGenerator(txPSDU,cfgHT); % generate txWaveform
  tx = [tx; zeros(15,cfgHT.NumTransmitAntennas)];  % Add trailing zeros to allow for channel filter delay
- **Transmitter-to-tag link**: The WiFi excitation signal propagates through the transmitter-to-tag channel to the tag. Since the transmitter and the tag are typically deployed in close proximity (on the order of several tens of centimeters), the transmitter-to-tag channel can be reasonably modeled as a flat-fading channel. This channel is represented by a complex coefficient $h=c+dj$, where $c$ and $d$ represent the real and imaginary components, respectively, and the magnitude is constrained by $0<|h|^{2}<1$.
  ```
  for chan_tx_tag_idx1 = 1:numTags
    bx_coeff_TxTag_real = (-1+(1+1)*rand(1,1)).*0.1;
    bx_coeff_TxTag_imag = (-1+(1+1)*rand(1,1)).*0.1;
    bx_coeff_TxTag = bx_coeff_TxTag_real + 1i*bx_coeff_TxTag_imag;
    tmp_exSig = tx.*bx_coeff_TxTag;
    exSig = [exSig,tmp_exSig];
  end
- **Tag data generation**:
  ```
  tagData = randi([0,1],numTagData-pLen,numTags);
  tagData = [ones(pLen,numTags);tagData];
- **Tag modulation**: The tag performs the backscatter operation. When bit 1 is transmitted, it introduces a phase of $\pi$; otherwise, no phase is introduced.
  ```
  for tag_idx1 = 1:numTags
    bxSig{tag_idx1} = survey_FreeCollision_funcBackscatter(exSig(:,tag_idx1),tagData(:,tag_idx1),1);
  end
- **Tag-to-receiver link**: The modulated signal then propagates through the tag-to-receiver channel to the receiver. The tag and the receiver are usually separated by a much larger distance, for which a multipath propagation model is more appropriate. Accordingly, we model the tag-to-receiver channel using MATLAB's built-in \emph{TGn Channel Model-B}, which captures the effects of indoor WiFi multipath propagation, CFO, and frequency selectivity.
  ```
  for chan_tag_rx_idx1 = 1:numTags
    reset(tgnChannel);
    rece_Sig{chan_tag_rx_idx1} = tgnChannel(bxSig{chan_tag_rx_idx1});
  end
- **Signal superposition from multiple tags**: The backscattered signal from multiple tags is superposed over the air. After that, noise is introduced.
  ```
  rxSig = complex(zeros(length(rece_Sig{1}),1));
  for rx_idx1 = 1:numTags
    rxSig = rxSig + rece_Sig{rx_idx1};
  end
  [rxSig,~,~] = func_awgn(rxSig,snr(i),'measured');
- **Demodulation**: The receiver decodes each tag’s data according to the corresponding system design.

# Packet Structure
The package consists of a funcs folder and executable MATLAB code files (in .m format). The functions in the funcs folder serve as function libraries required to support the execution of the .m files. All .m files can be executed directly in MATLAB R2021a.
- survey_ConcurScatter_Impact_of_PSDU_length.m file corresponds to Fig. 4 in the manuscript.
- survey_ConcurScatter_Impact_of_WiFi_transceiver_channel.m file corresponds to Fig. 6 in the manuscript.
- survey_ConcurScatter_Impact_of_backscatter_channel.m file corresponds to Fig. 7 in the manuscript.

# License

This project is for academic and research use only. Commercial use is prohibited.

# Modification

If you decide to use this code for academic issues, please support us by citing any of the following studies:

[1] Wu W, Xi W, Deng X, et al. Leveraging time-shifted orthogonal codes for concurrent backscatter communication[J]. IEEE Transactions on Sustainable Computing, 2024. DOI: 10.1109/TSUSC.2024.3517612. URL:https://ieeexplore.ieee.org/abstract/document/10803113

@article{wu2024leveraging,
  title={Leveraging time-shifted orthogonal codes for concurrent backscatter communication},
  author={Wu, Weiqi and Xi, Wei and Deng, Xianjun and Wang, Shuai and Zhou, Haoquan and Gong, Wei},
  journal={IEEE Transactions on Sustainable Computing},
  volume={10},
  number={4},
  pages={666--677},
  year={2024},
  publisher={IEEE}
}

[2] Wu W, Wang X, Hawbani A, et al. A survey on ambient backscatter communications: Principles, systems, applications, and challenges[J]. Computer Networks, 2022, 216: 109235.

@article{wu2022survey,
  title={A survey on ambient backscatter communications: Principles, systems, applications, and challenges},
  author={Wu, Weiqi and Wang, Xingfu and Hawbani, Ammar and Yuan, Longzhi and Gong, Wei},
  journal={Computer Networks},
  volume={216},
  pages={109235},
  year={2022},
  publisher={Elsevier}
}

[3] Wu W, Gong W. Concurrent WiFi backscatter communication using a single receiver in IoT networks[J]. Computer Networks, 2025, 258: 111029.

@article{wu2025concurrent,
  title={Concurrent WiFi backscatter communication using a single receiver in IoT networks},
  author={Wu, Weiqi and Gong, Wei},
  journal={Computer Networks},
  volume={258},
  pages={111029},
  year={2025},
  publisher={Elsevier}
}

