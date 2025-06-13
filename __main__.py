import numpy as np
from potential import get_potential
from simulation import run_simulation

p_init = 10
sigma_init = 0.2
x0_init = -2.0

# ! === WYBÓR SYMULACJI === ! #
# Zmień tę jedną linię, aby uruchomić inną symulację!
# ? Dostępne opcje: 'brak', 'bariera A', 'bariera B', 'bariera waska', 'bariera nieskonczona', 'skok', 'spadek' , 'studnia', 'studnia nieskonczona', 'trololololo'

POTENTIAL_CHOICE = 'trololololo'

# ! === KONIEC WYBORU SYMULACJI === ! #
    
if __name__ == '__main__':
    E_kin = p_init**2 / (2 * 1)
    x_grid = np.linspace(-5, 5, 1000)
    V_potential = get_potential(POTENTIAL_CHOICE, x_grid, E_kin, p_init)
    run_simulation(
        V=V_potential, 
        potential_name=POTENTIAL_CHOICE, 
        p=p_init, 
        sigma=sigma_init, 
        x0=x0_init, 
        xmin=-5, 
        xmax=5, 
        N=1000
    )