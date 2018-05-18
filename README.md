# Filter Bank Multicarrier Modulation (FBMC)


This repository compares FBMC to OFDM based schemes. All Figure from R. Nissel, S. Schwarz, and M. Rupp, [“Filter bank multicarrier modulation schemes for future mobile communications”](https://publik.tuwien.ac.at/files/publik_260162.pdf), IEEE Journal on Selected Areas in Communications, 2017, 
can be reproduce. In particular this repository

* calculates the Power Spectral Density (PSD),
* simulates a MIMO transmission, 
* calculates the Signal-to-Interference Ratio (SIR) in doubly-selective channels, 
* calculates the time-frequency efficiency,
* simulates the throughput.  


## Requirements
We used Windows 7 (64bit) and Matlab R2013b/2016a, but newer versions (and some older) should also work. Note that Figure 12 requires the Matlab “Communications System Toolbox” for turbo coding.

## Reproducible Figures
The figure numbers are the same as in  [“Filter bank multicarrier modulation schemes for future mobile communications”](https://publik.tuwien.ac.at/files/publik_260162.pdf):

* **Figure  1**: 
Just an illustration.	

* **Figure  2**: 
Please run [`Figure_02_PowerSpectralDensity.m`](Figure_02_PowerSpectralDensity.m).

* **Figure  3**: 
Please run [`Figure_03_BERoverSNR_MIMO.m`](Figure_03_BERoverSNR_MIMO.m).

* **Figure  4**: 
Just an illustration.

* **Figure  5**: 
Please run [`OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m`](OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m). Note that this script uses pre calculated values, stored in [`OptimalSubcarrierSpacing/Results/`](OptimalSubcarrierSpacing/Results/). To generate those pre-calculated values, the script [`OptimalSubcarrierSpacing/Calculate_SIR_SubcarrierSpacing_Velocity.m`](OptimalSubcarrierSpacing/Calculate_SIR_SubcarrierSpacing_Velocity.m) needs to be executed.

* **Figure  6**: 
Please run [`OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m`](OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m), see comment of Figure 5

* **Figure  7**: 
Please run [`OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m`](OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m), see comment of Figure 5

* **Figure  8**: 
Please run [`OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m`](OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m), see comment of Figure 5

* **Figure  9**: 
Please run [`Figure_09_SpectralEfficiency.m`](Figure_09_SpectralEfficiency.m).

* **Figure  10**: 
Please run [`Figure_10_11_TwoSubcarrierSpacingsSameBand.m`](Figure_10_11_TwoSubcarrierSpacingsSameBand.m).

* **Figure  11**: 
Please run [`Figure_10_11_TwoSubcarrierSpacingsSameBand.m`](Figure_10_11_TwoSubcarrierSpacingsSameBand.m).

* **Figure  12**: 
Please run [`Figure_12_Throughput.m`](Figure_12_Throughput.m).

* **Figure  13**: 
Please run [`Figure_13_PowerSpectralDensityQuantization.m`](Figure_13_PowerSpectralDensityQuantization.m).


## Additional Explanations of FBMC

1. [`Explained_A_PrototypeFilters.m`](Explained_A_PrototypeFilters.m): Illustration of different prototype filters.
2. [`Explained_B_FBMC_OQAM.m`](Explained_B_FBMC_OQAM.m): Describes a back-to-back FBMC-OQAM transmission, based on Section III. In particular, the transmit matrix (18)-(22) as well as the IFFT approach (32) are implemented.
3. [`Explained_C_Coded_FBMC_OQAM.m`](Explained_C_Coded_FBMC_OQAM.m): Shows how to find the precoding matrix, see (26) and (27), allowing QAM transmissions in FBMC-OQAM at full rate. Furthermore, it illustrates the time/frequency spreading concept.
4. [`Explained_D_SIR_DoublySelectiveChannel.m`](Explained_D_SIR_DoublySelectiveChannel.m): Implements Equation (35)-(40). Furthermore, it compares the theoretical values to simulations.


## 5G New Radio
We also include a comparison to New Radio (NR) waveforms (WOLA, UFMC, f-OFDM):	

1. [`NR_5G_BER_DoublySelectiveChannel.m`](NR_5G_BER_DoublySelectiveChannel.m): Simulates the Bit Error Ratio (BER) in a doubly-selective channel. It compares FBMC-OQAM, CP-OFDM, WOLA, UFMC and f-OFDM.
	
2. [`NR_5G_SIR_TimeFrequencyOffset.m`](NR_5G_SIR_TimeFrequencyOffset.m): Calculates the Signal-to-Interference Ratio (SIR) in case of a time and a frequency offset for FBMC-OQAM, FBMC-QAM, CP-OFDM, WOLA, UFMC and f-OFDM.	



## Please Cite Our Paper

    @ARTICLE{Nissel2017,
		author  = {R. Nissel and S. Schwarz and M. Rupp},
		journal = {IEEE Journal on Selected Areas in Communications},
		title   = {Filter Bank Multicarrier Modulation Schemes for Future Mobile Communications},
		year 	= {2017},
		volume 	= {35},
		number 	= {8},
		pages 	= {1768-1782}, 
		doi 	= {10.1109/JSAC.2017.2710022},
		ISSN 	= {0733-8716},
		month 	= {Aug},
	}


## References
- R. Nissel, S. Schwarz, and M. Rupp, [“Filter bank multicarrier modulation schemes for future mobile communications”](https://publik.tuwien.ac.at/files/publik_260162.pdf) IEEE Journal on Selected Areas in Communications, vol. 35, no. 8, pp. 1768–1782, 2017.
- R. Nissel, [“Filter bank multicarrier modulation for future wireless systems”](http://publik.tuwien.ac.at/files/publik_265168.pdf), Dissertation, TU Wien, 2017.



