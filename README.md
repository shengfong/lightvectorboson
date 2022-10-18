# Light vector boson model file for ACROPOLIS

When using this code, please cite our work

- **BBN Photodisintegration Constraints on Gravitationally Produced Vector Bosons**\
  Chee Sheng Fong, Moinul Hossain Rahat, Shaikh Saad\
  https://arxiv.org/abs/2206.02802
  
  as well as the work by the creators of the original code [ACROPOLIS](https://github.com/hep-mh/acropolis).
  
## Abstract

We implement the model of a general light vector boson $V$ with mass in the range MeV to GeV for the public code ACROPOLIS 
to study the photodisintegration effects of light elements produced from the Big Bang Nucleosynthesis when the Universe
was about $10^4$ s old. The primary photon and electron/positron spectra from all five possible channels for $m_V \leq 1$ GeV are computed: $V \to e^+ e^-$, $V \to \mu^+ \mu^-$, $V \to \pi^+ \pi^-$, $V \to \pi^0 \gamma$, and $V \to \pi^0 \pi^+ \pi^-$. We have made appropriate modifications to ACROPOLIS to be able to take into account two monochromatic injection energies: that of electron/positron in $V \to e^+ e^-$ and that of photon in $V \to \pi^0 \gamma$. The model file  ``decay_vector_model.py`` contains the analytic expressions for the primary photon and electron/positron spectra for $V \to e^+e^-$ and $V \to \pi^0 \gamma$ as well as the final state radiation of photons for $V\to \mu^+\mu^-$ and $V\to \pi^+ \pi^-$. The primary photon and electron/positron spectra for $V\to \mu^+\mu^-$, $V\to \pi^+ \pi^-$, and $V \to \pi^0 \pi^+ \pi^-$ are precalculated numerically and the data is stored in the folder ``\spec_data`` and will be called automatically in the calculations. By default, the Boltzmann equations are only solved up till the epoch of matter-radiation equality when the cosmic time is $t = 2\times 10^{12}$ s. This means that for a very long-lived vector boson particle $\tau \gtrsim 10^{12}$ s, the bounds obtained are conservative. 

## Changelog

**June 4, 2022** \
From the version 1.2.1 of [ACROPOLIS](https://github.com/hep-mh/acropolis), changes were made to ``cascade.py``, ``models.py``, ``nucl.py`` and ``params.py`` to adapt to our model file ``decay_vector_model.py``. For more details about this model, please refer to our [work](https://arxiv.org/abs/22xx.xxxxx) and for more detail on the original code, please refer to [ACROPOLIS](https://github.com/hep-mh/acropolis).



## Usage

To verify if the parameter space of a light vector boson is viable, one can run the executable ``decayvector`` in the terminal with the command

```
./decayvector 700 1e8 1 1e-6 0.1 0.1 0.7 0.01 0.09
```
where the input parameters after ``./decayvector`` are

- $m_V$ [MeV]: light vector mass 
- $\tau$ [s]: light vector lifetime
- $T_0$ [MeV]: reference cosmic temperature
- $\frac{n_V}{n_\gamma}(T_0)$: the number density of $V$ over the photon density $n_V/n_\gamma$ at $T_0$
- $BR_{ee}$: decay branching ratio of $V \to e^+ e^-$
- $BR_{\mu\mu}$: decay branching ratio of $V \to \mu^+ \mu^-$
- $BR_{\pi\pi}$: decay branching ratio of $V \to \pi^+ \pi^-$
- $BR_{\pi\gamma}$: decay branching ratio of $V \to \pi^0 \gamma$
- $BR_{3\pi}$: decay branching ratio of $V \to \pi^0 \pi^+ \pi^-$

The number density of $V$ over the photon density $n_V/n_\gamma$ is related to $Y_V = n_V/s$ where $s$ is the cosmic entropy density as follows

$$\frac{n_V}{n_\gamma}(T_0) = \frac{\pi^4 g_\star(T_0)}{45\zeta(3)} Y_V.$$

By default, the combined exclusion at 95\% CL for Helium-4 and deuterium considering the sum of theoretical and experimental errors in quadrature (see https://arxiv.org/abs/2206.02802 for further details). An example of the output is as follows 
```
Results:   Yp = 0.224415, H2/p = 0.000395, He3/p = 0.006490
Excluded by the BBN measurements at 2 sigma 
(default: He3/p not considered).
Runtime - - - 56.018600 mins - - -
```
``Yp``, ``H2/p``, and ``He3/p`` denotes respectively the Helium-4 mass fraction, $n_{\rm D}/n_{\rm H}$, and $n_{\rm He^3}/n_{\rm H}$.

We also include ``decayvectorscan`` to illustrate how to scan the parameter space in mass $m_V$ and lifetime $\tau$ of a light vector boson. The branching ratios for the decay of a dark photon which only couples to the Standard Model photon through kinetic mixing is used in this example (the branching ratios are read from ``spec_data/darkphoton_MV_BR.dat``).

