===========================================
If you use this code, please cite our paper
===========================================
@ARTICLE{Nissel2017,
	author  = {R. Nissel and S. Schwarz and M. Rupp},
	journal = {IEEE Journal on Selected Areas in Communications},
	title   = {Filter Bank Multicarrier Modulation Schemes for Future Mobile Communications},
	year	= {2017, to appear},
	month 	= {},
}


=====================
System Requirements
===================== 
We used Windows 7 (64bit) and MATLAB R2013b/2016a, but newer versions (and some older) should also work.

Figure 12 requires the "Communications System Toolbox" for the turbo coding.



====================
Description
====================

All of our figures can be reproduced:
	Figure  1:	  Just an illustration.	
	Figure  2: 	Please run "Figure_02_PowerSpectralDensity.m".
	Figure  3: 	Please run "Figure_03_BERoverSNR_MIMO.m".
	Figure  4: 	  Just an illustration.
	Figure  5:	Please run "OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m". Note that this script uses pre-calculated values which are stored in "OptimalSubcarrierSpacing/Results/". To generate those pre-calculated values, the script "OptimalSubcarrierSpacing/Calculate_SIR_SubcarrierSpacing_Velocity.m" needs to be executed.
	Figure  6: 	Please run "OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m", see above.
	Figure  7: 	Please run "OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m", see above.
	Figure  8: 	Please run "OptimalSubcarrierSpacing/Figure_5_6_7_8_SIR_OptimalSubcarrierSpacing.m", see above.
	Figure  9: 	Please run "Figure_09_SpectralEfficiency.m".
	Figure 10:	Please run "Figure_10_11_TwoSubcarrierSpacingsSameBand.m".
	Figure 11:	Please run "Figure_10_11_TwoSubcarrierSpacingsSameBand.m".
	Figure 12: 	Please run "Figure_12_Throughput.m".
	Figure 13:	Please run "Figure_13_PowerSpectralDensityQuantization.m".


We also include additional explanations of FBMC
	1) "Explained_A_PrototypeFilters.m"			This script plots different prototype filters.
	2) "Explained_B_FBMC_OQAM.m"				This script describes a back-to-back FBMC-OQAM transmission, based on Section III. In particular, the transmit matrix (18)-(22) as well as the IFFT approach (32) are implemented.
	3) "Explained_C_Coded_FBMC_OQAM.m"			This script shows how to find the precoding matrix, see (26) and (27), which allows QAM transmission in FBMC-OQAM at full rate. Furthermore, it illustrates the time/frequency spreading concept.
	4) "Explained_D_SIR_DoublySelectiveChannel.m"		This script implements Equation (35)-(40). Furthermore, it compares the theoretical values to simulations.


Furthermore, we include a comparison of New Ratio (NR) waveforms (WOLA, UFMC, f-OFDM):	
	1) "NR_5G_BER_DoublySelectiveChannel.m"			This script simulates the Bit Error Ratio (BER) in a doubly-selective channel. It compares FBMC-OQAM, CP-OFDM, WOLA, UFMC, f-OFDM.
	2) "NR_5G_SIR_TimeFrequencyOffset.m"			This script calculates the Signal-to-Interference Ratio (SIR) in case of a time	and a frequency offset for FBMC-OQAM, FBMC-QAM, CP-OFDM, WOLA, UFMC and f-OFDM.	




