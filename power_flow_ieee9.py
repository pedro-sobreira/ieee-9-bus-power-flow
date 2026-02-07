import numpy as np
import os

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

def power_flow():
    clear_screen()
    print("====================================================")
    print("Cálculo de Fluxo de Potência - IEEE 9 Bus System")
    print("Método: Newton-Raphson")
    print("====================================================\n")

    # 1. Configurações do Sistema
    base_mva = 100.0
    max_iter = 20
    tolerance = 1e-6

    # 2. Dados dos Barramentos (Bus Data)
    # Tipo: 1=Slack, 2=PV, 3=PQ
    # [Bus, Type, V_pu, Angle_deg, P_gen_mw, Q_gen_mvar, P_load_mw, Q_load_mvar]
    bus_data = np.array([
        [1, 1, 1.040, 0.0, 0.0, 0.0, 0.0, 0.0],
        [2, 2, 1.025, 0.0, 163.0, 0.0, 0.0, 0.0],
        [3, 2, 1.025, 0.0, 85.0, 0.0, 0.0, 0.0],
        [4, 3, 1.000, 0.0, 0.0, 0.0, 0.0, 0.0],
        [5, 3, 1.000, 0.0, 0.0, 0.0, 125.0, 50.0],
        [6, 3, 1.000, 0.0, 0.0, 0.0, 90.0, 30.0],
        [7, 3, 1.000, 0.0, 0.0, 0.0, 0.0, 0.0],
        [8, 3, 1.000, 0.0, 0.0, 0.0, 100.0, 35.0],
        [9, 3, 1.000, 0.0, 0.0, 0.0, 0.0, 0.0]
    ])

    # 3. Dados das Linhas e Transformadores (Line Data)
    # [From, To, R_pu, X_pu, B_half_pu]
    line_data = np.array([
        [1, 4, 0.0000, 0.0576, 0.0000],
        [4, 5, 0.0100, 0.0850, 0.0880],
        [4, 6, 0.0170, 0.0920, 0.0790],
        [5, 7, 0.0320, 0.1610, 0.1530],
        [2, 7, 0.0000, 0.0625, 0.0000],
        [7, 8, 0.0085, 0.0720, 0.0745],
        [8, 9, 0.0119, 0.1008, 0.1045],
        [3, 9, 0.0000, 0.0586, 0.0000],
        [6, 9, 0.0390, 0.1700, 0.1790]
    ])

    num_buses = len(bus_data)
    
    # 4. Construção da Matriz Ybus
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    for line in line_data:
        f = int(line[0]) - 1
        t = int(line[1]) - 1
        r = line[2]
        x = line[3]
        b_half = line[4]
        
        z = r + 1j*x
        y = 1/z
        
        Ybus[f, f] += y + 1j*b_half
        Ybus[t, t] += y + 1j*b_half
        Ybus[f, t] -= y
        Ybus[t, f] -= y

    # 5. Inicialização de Variáveis
    V = bus_data[:, 2].copy()
    theta = np.radians(bus_data[:, 3].copy())
    P_spec = (bus_data[:, 4] - bus_data[:, 6]) / base_mva
    Q_spec = (bus_data[:, 5] - bus_data[:, 7]) / base_mva
    
    # Identificar tipos de barramento
    slack_idx = np.where(bus_data[:, 1] == 1)[0]
    pv_idx = np.where(bus_data[:, 1] == 2)[0]
    pq_idx = np.where(bus_data[:, 1] == 3)[0]
    
    non_slack_idx = np.concatenate([pv_idx, pq_idx])
    
    # 6. Iterações de Newton-Raphson
    for iteration in range(max_iter):
        # Calcular P e Q calculados
        P_calc = np.zeros(num_buses)
        Q_calc = np.zeros(num_buses)
        
        for i in range(num_buses):
            for j in range(num_buses):
                P_calc[i] += V[i] * V[j] * (Ybus[i, j].real * np.cos(theta[i] - theta[j]) + 
                                           Ybus[i, j].imag * np.sin(theta[i] - theta[j]))
                Q_calc[i] += V[i] * V[j] * (Ybus[i, j].real * np.sin(theta[i] - theta[j]) - 
                                           Ybus[i, j].imag * np.cos(theta[i] - theta[j]))
        
        # Calcular resíduos
        dP = P_spec[non_slack_idx] - P_calc[non_slack_idx]
        dQ = Q_spec[pq_idx] - Q_calc[pq_idx]
        
        mismatch = np.concatenate([dP, dQ])
        if np.max(np.abs(mismatch)) < tolerance:
            print(f"Convergência atingida em {iteration} iterações.\n")
            break
            
        # 7. Construção da Matriz Jacobiana
        # J = [[H, N], [M, L]]
        # H = dP/dtheta, N = dP/dV, M = dQ/dtheta, L = dQ/dV
        
        # Índices para a Jacobiana
        n_ns = len(non_slack_idx)
        n_pq = len(pq_idx)
        J = np.zeros((n_ns + n_pq, n_ns + n_pq))
        
        # Preencher H (dP/dtheta)
        for i_idx, i in enumerate(non_slack_idx):
            for j_idx, j in enumerate(non_slack_idx):
                if i == j:
                    J[i_idx, j_idx] = -Q_calc[i] - (V[i]**2) * Ybus[i, i].imag
                else:
                    J[i_idx, j_idx] = V[i] * V[j] * (Ybus[i, j].real * np.sin(theta[i] - theta[j]) - 
                                                    Ybus[i, j].imag * np.cos(theta[i] - theta[j]))
        
        # Preencher N (dP/dV)
        for i_idx, i in enumerate(non_slack_idx):
            for j_idx, j in enumerate(pq_idx):
                if i == j:
                    J[i_idx, n_ns + j_idx] = (P_calc[i] / V[i]) + Ybus[i, i].real * V[i]
                else:
                    J[i_idx, n_ns + j_idx] = V[i] * (Ybus[i, j].real * np.cos(theta[i] - theta[j]) + 
                                                    Ybus[i, j].imag * np.sin(theta[i] - theta[j]))
        
        # Preencher M (dQ/dtheta)
        for i_idx, i in enumerate(pq_idx):
            for j_idx, j in enumerate(non_slack_idx):
                if i == j:
                    J[n_ns + i_idx, j_idx] = P_calc[i] - (V[i]**2) * Ybus[i, i].real
                else:
                    J[n_ns + i_idx, j_idx] = -V[i] * V[j] * (Ybus[i, j].real * np.cos(theta[i] - theta[j]) + 
                                                            Ybus[i, j].imag * np.sin(theta[i] - theta[j]))
        
        # Preencher L (dQ/dV)
        for i_idx, i in enumerate(pq_idx):
            for j_idx, j in enumerate(pq_idx):
                if i == j:
                    J[n_ns + i_idx, n_ns + j_idx] = (Q_calc[i] / V[i]) - Ybus[i, i].imag * V[i]
                else:
                    J[n_ns + i_idx, n_ns + j_idx] = V[i] * (Ybus[i, j].real * np.sin(theta[i] - theta[j]) - 
                                                           Ybus[i, j].imag * np.cos(theta[i] - theta[j]))
        
        # Resolver sistema linear
        dx = np.linalg.solve(J, mismatch)
        
        # Atualizar variáveis
        theta[non_slack_idx] += dx[:n_ns]
        V[pq_idx] += dx[n_ns:]
    else:
        print("O método não convergiu no número máximo de iterações.")

    # 8. Resultados Finais
    print("Resultados do Fluxo de Potência:")
    print("-" * 65)
    print(f"{'Bus':<5} {'V (pu)':<10} {'Ângulo (deg)':<15} {'P (MW)':<12} {'Q (MVAr)':<12}")
    print("-" * 65)
    
    # Recalcular P e Q para todos os barramentos para exibição
    P_final = np.zeros(num_buses)
    Q_final = np.zeros(num_buses)
    for i in range(num_buses):
        for j in range(num_buses):
            P_final[i] += V[i] * V[j] * (Ybus[i, j].real * np.cos(theta[i] - theta[j]) + 
                                       Ybus[i, j].imag * np.sin(theta[i] - theta[j]))
            Q_final[i] += V[i] * V[j] * (Ybus[i, j].real * np.sin(theta[i] - theta[j]) - 
                                       Ybus[i, j].imag * np.cos(theta[i] - theta[j]))
        
        print(f"{i+1:<5} {V[i]:<10.4f} {np.degrees(theta[i]):<15.4f} {P_final[i]*base_mva:<12.4f} {Q_final[i]*base_mva:<12.4f}")
    print("-" * 65)

if __name__ == "__main__":
    power_flow()
