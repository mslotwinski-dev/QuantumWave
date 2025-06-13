import numpy as np

def get_potential(name, x, E_kinetyczna, p):
    V = np.zeros_like(x)

    print(f"Wybrano potencjał: '{name}'")

    if name == 'brak': 
        pass

    elif name == 'bariera A':
        V[(x > 0)] = E_kinetyczna

    elif name == 'bariera B':
        V[(x < 0)] = E_kinetyczna

    elif name == 'bariera waska':
        V[(x > 0) & (x < 1)] = E_kinetyczna * 10000

    elif name == 'bariera nieskonczona':
        V[(x > 0) & (x < 1)] = E_kinetyczna * 10000

    elif name == 'skok':
        V[(x < 0)] = 0.75*E_kinetyczna
        V[(x > 0)] = 1.5*E_kinetyczna

    elif name == 'spadek':
        V = E_kinetyczna *(1- 0.2*x)
  
    elif name == 'studnia':
        V[(x < -3) | (x > -1)] = E_kinetyczna

    elif name == 'studnia nieskonczona':
        V[(x < -3) | (x > -1)] = E_kinetyczna * 10000

    elif name == 'trololololo':
        V = E_kinetyczna * np.sin(x)**2 - 0.5 * E_kinetyczna + E_kinetyczna *(1- 0.2*x**2)
        
    else:
        raise ValueError(f"Nieznany typ potencjału: {name}. Dostępne opcje: 'brak', 'bariera A', 'bariera B', 'bariera waska', 'bariera nieskonczona', 'skok', 'spadek' , 'studnia', 'studnia nieskonczona', 'trololololo'")
        
    return V
