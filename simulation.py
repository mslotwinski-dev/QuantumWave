import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def run_simulation(V, potential_name, p, sigma, x0, m=1, hbar=1, xmin=-5, xmax=5, N=1000):
    # --- Konfiguracja siatki ---
    x = np.linspace(xmin, xmax, N)
    dx = x[1] - x[0]
    E_kinetyczna = p**2 / (2 * m)

    # --- Początkowa funkcja falowa ---
    x_wew = x[1:-1]
    Psi0 = np.exp(-(x_wew - x0)**2 / (4 * sigma**2)) * np.exp(1j * p * x_wew)
    A = np.sum(np.abs(Psi0)**2 * dx)
    Psi0 = Psi0 / np.sqrt(A)

    # --- Hamiltonian i jego diagonalizacja ---
    kinetic = (-hbar**2 / (2 * m * dx**2)) * (
        np.diag(np.ones(N - 3), -1) - 2 * np.diag(np.ones(N - 2)) + np.diag(np.ones(N - 3), 1)
    )
    potential = np.diag(V[1:-1])
    H = kinetic + potential

    print("Diagonalizuję Hamiltonian... (to może chwilę potrwać)")
    E, psi_T = np.linalg.eigh(H)
    psi = psi_T.T
    print("Diagonalizacja zakończona.")
    
    # --- Normalizacja funkcji falowych ---
    A_psi = np.sum(np.abs(psi[0])**2 * dx)
    psi = psi / np.sqrt(A_psi)

    # --- Obliczanie współczynników rozwinięcia ---
    c = np.zeros(len(E), dtype=np.complex128)
    for i in range(len(E)):
        c[i] = np.sum(np.conj(psi[i]) * Psi0 * dx)

    # --- Konfiguracja wykresu i animacji ---
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.set_xlabel("x")
    ax.set_ylabel("Amplituda / Gęstość prawdopodobieństwa")
    ax.set_title(f"Ewolucja paczki falowej - Potencjał: {potential_name.replace('_', ' ').title()}")
    ax.grid(True)

    max_initial_prob_density = np.max(np.abs(Psi0)**2)
    V_max = np.max(np.abs(V))
    if V_max > 0:
        scaling_factor = max_initial_prob_density / V_max
        ax.plot(x, V * scaling_factor, color='red', label='Potencjał V(x) (przeskalowany)', lw=1.5, alpha=0.7)

    line_prob, = ax.plot([], [], color='#00647D', lw=2.5, label='$|\Psi(x,t)|^2$')
    line_real, = ax.plot([], [], color='#CE9F6D', lw=1.5, alpha=0.8, label='Re($\Psi$)')
    line_imag, = ax.plot([], [], color='#CE9FD5', lw=1.5, alpha=0.8, label='Im($\Psi$)')
    ax.legend()
    
    max_y = np.max(np.abs(Psi0))
    ax.set_ylim(-max_y * 1.2, max_initial_prob_density * 1.2)
    ax.set_xlim(xmin, xmax)

    def init():
        line_prob.set_data([], [])
        line_real.set_data([], [])
        line_imag.set_data([], [])
        return line_prob, line_real, line_imag,

    def update(frame):
        t = frame * 0.001
        Psi_wew = np.zeros(N - 2, dtype=np.complex128)
        for i in range(len(E)):
            Psi_wew += c[i] * psi[i] * np.exp(-1j * E[i] * t / hbar)

        full_prob_density = np.zeros(N)
        full_real_part = np.zeros(N)
        full_imag_part = np.zeros(N)
        
        full_prob_density[1:-1] = np.abs(Psi_wew)**2
        full_real_part[1:-1] = np.real(Psi_wew)
        full_imag_part[1:-1] = np.imag(Psi_wew)
        
        line_prob.set_data(x, full_prob_density)
        line_real.set_data(x, full_real_part)
        line_imag.set_data(x, full_imag_part)
        
        return line_prob, line_real, line_imag,

    ani = FuncAnimation(fig, update, frames=400, init_func=init, blit=True, interval=20)
    plt.show()
