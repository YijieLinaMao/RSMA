# Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA.

This is a code package related to the following paper, which receives the Best Paper Award of [EURASIP JWCN 2022](https://eurasip.org/newsletter/newsletter_2022-04.html):

Y. Mao, B. Clerckx and V. O. K. Li, "[Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA](https://link.springer.com/article/10.1186/s13638-018-1104-7)," EURASIP Journal on Wireless Communications and Networking 2018.1 (2018): 133.

# Content of Code Package
Here is a detailed description of the package: 
- The code in all packages are implemented in Matlab environment with CVX toolbox assisted. 
- In the 'Figure 6b' code package,  Fig. 6(b) of the above paper will be reproduced by running the Matlab script 'main.m'. By changing the variables 'bias' (channel gain difference between the users), 'NT'(number of transmit antenna), you can reproduce Fig. 5--Fig. 6.
- In the 'Figure 7' code package,  Fig. 7 of the above paper will be reproduced by running the Matlab script 'main.m'. By changing the variable 'bias' and the channel realizations, you can reproduce Fig. 7--Fig. 10.
- In the 'Figure 11' code package,  Fig. 11 of the above paper will be reproduced by running the Matlab script 'main.m'. By changing the variable 'bias' and the channel realizations, you can reproduce Fig. 11--Fig. 12.
- In the 'Figure 13' code package, Fig. 13 of the above paper will be reproduced by running the Matlab script 'main.m'. By changing the variable 'weight' and 'NT', you can reproduce Fig. 13, Fig. 14, Fig. 16.

# Abstract of the Article
Space-division multiple access (SDMA) utilizes linear precoding to separate users in the spatial domain and relies on fully treating any residual multi-user interference as noise. Non-orthogonal multiple access (NOMA) uses linearly precoded superposition coding with successive interference cancellation (SIC) to superpose users in the power domain and relies on user grouping and ordering to enforce some users to fully decode and cancel interference created by other users.

In this paper, we argue that to efficiently cope with the high throughput, heterogeneity of quality of service (QoS), and massive connectivity requirements of future multi-antenna wireless networks, multiple access design needs to depart from those two extreme interference management strategies, namely fully treat interference as noise (as in SDMA) and fully decode interference (as in NOMA).

Considering a multiple-input single-output broadcast channel, we develop a novel multiple access framework, called rate-splitting multiple access (RSMA). RSMA is a more general and more powerful multiple access for downlink multi-antenna systems that contains SDMA and NOMA as special cases. RSMA relies on linearly precoded rate-splitting with SIC to decode part of the interference and treat the remaining part of the interference as noise. This capability of RSMA to partially decode interference and partially treat interference as noise enables to softly bridge the two extremes of fully decoding interference and treating interference as noise and provides room for rate and QoS enhancements and complexity reduction.

The three multiple access schemes are compared, and extensive numerical results show that RSMA provides a smooth transition between SDMA and NOMA and outperforms them both in a wide range of network loads (underloaded and overloaded regimes) and user deployments (with a diversity of channel directions, channel strengths, and qualities of channel state information at the transmitter). Moreover, RSMA provides rate and QoS enhancements over NOMA at a lower computational complexity for the transmit scheduler and the receivers (number of SIC layers).


# License and Referencing
This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

# Acknowledgements
The paper was supported by the UK Engineering and Physical Sciences Research Council (EPSRC) under grant EP/N015312/1.

