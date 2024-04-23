# Theory

## Single-qubit gates
On this page we are describing how a single-qubit gate that implements rotation around the x- or y-axis on the Bloch sphere can be performed.

We consider a driven weakly anharmonic qubit whose Hamiltonian in lab frame can be written as

\begin{equation}
    \label{eq:Transmon}
    \frac{H}{\hbar} = \omega_q a^\dagger a+\frac{\alpha}{2} a^\dagger a^\dagger a a +\mathcal{E}(t)a^\dagger+\mathcal{E}(t)^*a,
\end{equation}

where $\omega_q\equiv \omega_q^{0\rightarrow 1}$ is the qubit frequency and $\alpha = \omega_q ^{1\rightarrow 2}-\omega_q^{0\rightarrow 1}$ is the anharmonicity. The driving and control is given by

\begin{equation}
    E(t)= \begin{cases} 
        \Omega^x(t)\cos(\omega_d t)+\Omega^y(t)\sin(\omega_d t),& 0<t<t_g, \\ 0, & \text{otherwise}.
    \end{cases}
\end{equation}

Here $\Omega^x(t)$ and $\Omega^y(t)$ are two independent quadrature controls, $t_g$ is the total gate-time, and $\omega_d$ is the drive frequency. Next we move into the rotating frame of the drive by performing the following unitary transformation $U(t)=e^{i\omega_r t a^\dagger a}$, where $\omega_r$ is the rotating frame frequency. The Hamiltonian in the rotating frame after having performed the rotating wave approximation reads

\begin{multline}
    \frac{H^R}{\hbar}=
    \Delta a^\dagger a + \frac{\alpha}{2} a^{\dagger 2}a^2 + 
    (\frac{\Omega^x(t)}{2}\cos([\omega_r-\omega_d]t)-\frac{\Omega^y(t)}{2}\sin([\omega_r-\omega_d]t))(a^\dagger + a) \\
    + (\frac{\Omega^x(t)}{2}\sin([\omega_r-\omega_d]t)+\frac{\Omega^y(t)}{2}\cos([\omega_r-\omega_d]t))(ia^\dagger - ia),
\end{multline}

where $\Delta \equiv \omega_q - \omega_r$ is the qubit detuning. 

As a concrete example, assume that we apply a pulse at the qubit frequency $\omega_d=\omega_q$, and choose the rotating frame of the drive $\omega_r=\omega_d$. Then,

\begin{equation}
    \frac{H^R}{\hbar} =
    \frac{\alpha}{2} a^{\dagger 2}a^2
    + \frac{\Omega^x(t)}{2}(a^\dagger + a)
    + \frac{\Omega^y(t)}{2}(ia^\dagger - ia).
\end{equation}

If we treat the Hamiltonian as an effective two level system (ignoring the anharmonic term) and make the replacement $(a^\dagger + a)\rightarrow \sigma_x$ and $(ia^\dagger-ia)\rightarrow \sigma_y$, we obtain

\begin{equation}
    \frac{H^R}{\hbar} = \frac{\Omega^x(t)}{2}\sigma_x + \frac{\Omega^y(t)}{2}\sigma_y,
\end{equation}

showing that an in-phase pulse (i.e. the $\Omega^x(t)$ quadrature component) corresponds to a rotation around the $x$-axis while the out-of-phase pulse (i.e. the $\Omega^y(t)$ quadrature component), corresponds to rotations about the $y$-axis. As a concrete example of an in-phase pulse, writing out the unitary evolution operator yields,

\begin{equation}
    U^R(t)=\exp([-\frac{i}{2}\int_0^t\Omega^x(t')\mathrm{d}t']\sigma_x).
\end{equation}

By defining the angle

\begin{equation}
    \Theta(t)=\int_0^t\Omega^x(t')\mathrm{d}t',
\end{equation}

which is the angle a state is rotated given a waveform envelope $\Omega^x(t)$. This means that to implement a $\pi$-pulse on the $x$-axis one would solve $\Theta(t)=\pi$ and output the signal in-phase with the qubit drive.

In this simple example we assumed that we could ignore the higher levels of the qubit. In general leakage errors which take the qubit out of the computational subspace as well as phase errors can occur. To combat theses errors the so-called DRAG[@motzoi2009simple] procedure (Derivative Reduction by Adiabatic Gate) is used. In doing so we apply an extra signal in the out-of-phase component, such that

\begin{align}
    \Omega^x(t) = B e^{-\frac{(t-t_g)^2}{2\sigma^2}},\quad
    \Omega^y(t) = q\sigma\frac{d\Omega^x(t)}{dt}
\end{align}

where $q$ is a scale parameter that needs to be optimized with respect to a $\pi/2$-pulse. Interchanging $\Omega^x(t)$ and $\Omega^y(t)$ in the equation above corresponds to DRAG pulsing the $\Omega^y(t)$ component. The amplitude $B$ is fixed such that

\begin{equation}
    \Big|\int_{0}^{t}[\Omega^x(t')+i\Omega^y(t')]\mathrm{d}t'\Big|=\pi.
\end{equation}

for a $\pi$-pulse with DRAG.