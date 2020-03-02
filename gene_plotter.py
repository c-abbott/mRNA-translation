import numpy as np 
import matplotlib.pyplot as plt

"""
    Dictionary containing data collected from
    early_time.py for a set of unique E. coli
    genes.
"""
gene_dict = {
            'yaaA': [15.1658, 13.2024, 228], 'talB': [15.1884, 13.4269, 204],
            'dnaK': [32.8133, 28.4524, 144], 'dnaJ': [20.9123, 17.8550, 150],
            'nhaA': [23.9344, 18.0158, 396], 'araC': [19.5319, 15.8746, 294],
            'araE': [28.0273, 22.7060, 498], 'asnB': [30.2830, 27.2119, 192],
            'mutS': [46.2525, 39.5050, 330], 'mglB': [18.7629, 15.8855, 624],
            'yeiE': [16.5390, 13.7827, 324], 'yeiH': [20.7176, 16.8940, 396],
            'yeiP': [11.0230, 9.0845, 228],  'menA': [18.1965, 14.8095, 1146],
            'ygiW': [7.3537, 6.0203, 954],   'ygiM': [12.3366, 10.7214, 228],
            'dapB': [15.6233, 13.0669, 210], 'leuL': [1.5813, 1.3545, 420],
            'secA': [45.7775, 41.9156, 156], 'yhaH': [7.3976, 6.0929, 384],
            'rpsL': [6.3650, 5.3821, 480],   'yacC': [7.3487, 6.0256, 414],
            'queA': [20.6838, 17.3114, 954], 'phoB': [13.7744, 11.2784, 294],
            'ahpC': [9.1262, 8.7535, 1086],  'cstA': [41.3607, 33.8641, 276],
            'nagE': [36.8743, 30.9129, 552], 'glnS': [30.8389, 26.0784, 144]
            }

# Repackage for plotting.
gene_data = {"T1_exp": [], "T1_ana": [], "t_1/2": [], "label": []}
for label, coord in gene_dict.items():
    gene_data["T1_exp"].append(coord[0])
    gene_data["T1_ana"].append(coord[1])
    gene_data["t_1/2"].append(coord[2])
    gene_data["label"].append(label)

# Statistiacal manipulation.
half_lives = np.array(gene_data["t_1/2"])
T1s_exp = np.array(gene_data["T1_exp"])
T1s_ana = np.array(gene_data["T1_ana"])
M = np.vstack([half_lives, np.ones(len(half_lives))]).T
m_exp, c_exp = np.linalg.lstsq(M, T1s_exp, rcond=None)[0]
m_ana, c_ana = np.linalg.lstsq(M, T1s_ana, rcond=None)[0]
pear_R_exp = np.corrcoef(half_lives, T1s_exp)[1, 0]
pear_R_ana = np.corrcoef(half_lives, T1s_exp)[1, 0]

# Display scatter plot data.
plt.figure(figsize=(10, 8))
plt.title('Translation Time vs. mRNA Half-Life', fontsize=20)
plt.xlabel('mRNA Half-Life (t_1/2)', fontsize=15)
plt.ylabel('Mean Translation Time (<T1>)', fontsize=15)
plt.scatter(gene_data["t_1/2"], gene_data["T1_exp"], marker='o', label = "Simulation Data")
plt.scatter(gene_data["t_1/2"], gene_data["T1_ana"], marker='1', label = "Analytical Data")
plt.plot(half_lives, m_exp*half_lives + c_exp, 'r',
         label="Simulation Fit |%6s, r = %6.2e" % ('red', pear_R_exp))
plt.plot(half_lives, m_ana*half_lives + c_ana, 'g',
         label="Analytical Fit |%6s, r = %6.2e" % ('green', pear_R_ana))
plt.legend()
plt.savefig("gene_plot.png")
plt.show()
