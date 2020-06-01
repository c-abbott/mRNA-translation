import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.style as style
import csv

"""
    Dictionary containing data collected from
    early_time.py for a set of unique E. coli
    genes.

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
            'nagE': [36.8743, 30.9129, 552], 'glnS': [30.8389, 26.0784, 144],
            'fldA': [9.7444, 8.1320, 162],   'ybfE': [5.1051, 4.1962, 264],
            'ybfF': [14.9854, 12.8428, 180], 'ybfA': [3.9918, 3.2873, 384],
            'gltA': [22.7271, 21.4408, 138], 'aroG': [20.0883, 17.0728, 156],
            'gpmA': [11.8463, 10.4654, 174], 'ybhC': [24.8642, 20.6214, 384],
            'ybhB': [9.9909, 8.4315, 324],   'bioA': [27.1148, 22.7728, 132],
            'ybiC': [22.8083, 19.9431, 228], 'ybiJ': [4.6561, 4.1471, 228],
            'glnH': [13.7425, 11.7472, 342], 'ybiT': [28.9572, 24.9153, 306],
            'ybjG': [12.3880, 10.5454, 492], 'rimK': [17.2669, 14.0255, 264],
            'ybjN': [9.9647, 8.2829, 324],   'potF': [21.8504, 18.5763, 252],
            'artJ': [14.4017, 13.0224, 246], 'cspD': [4.8501, 3.9999, 924],
            'infA': [3.8718, 3.2838, 204],   'trxB': [17.9053, 14.8329, 234],
            'serS': [23.3094, 19.4853, 156], 'ycaD': [23.2843, 19.0645, 324],
            'pflA': [14.4859, 12.2681, 258], 'pflB': [39.7954, 35.8073, 348],
            'serC': [20.6774, 17.7228, 120], 'aroA': [22.4971, 19.6186, 114],
            'rpsA': [26.6926, 22.4792, 252], 'aspC': [20.7294, 19.0304, 258]
            }
"""


def main():
    # Open file.
    with open("overall_gene_dat.csv", "r") as f:
        reader = csv.reader(f)
        # Skip first line.
        next(f)
        # Dictionary for plotting.
        gene_dict = {}
        # Assigning values from CSV.
        for row in reader:
            key = row[0]
            # Duplicate row handling.
            if key in gene_dict:
                pass
            # Making data floats.
            vals = row[1:]
            gene_dict[key] = [float(i) for i in vals]

    # Repackage for plotting.
    gene_data = {"T1_exp": [], "T1_ana": [],
                 "lambda": [], "t_1/2": [], "total_time": [], "label": []}
    for label, coord in gene_dict.items():
        gene_data["T1_exp"].append(coord[0])
        gene_data["T1_ana"].append(coord[1])
        gene_data["lambda"].append(coord[2])
        gene_data["t_1/2"].append(coord[3])
        gene_data["total_time"].append(coord[4])
        gene_data["label"].append(label)

    # Statistiacal manipulation.
    half_lives = np.array(gene_data["t_1/2"])
    T1s_exp = np.array(gene_data["total_time"])
    lambdas = np.array(gene_data["lambda"])
    M = np.vstack([half_lives, np.ones(len(half_lives))]).T
    m_exp, c_exp = np.linalg.lstsq(M, T1s_exp, rcond=None)[0]
    m_l, c_l = np.linalg.lstsq(M, lambdas, rcond=None)[0]
    pear_R_exp, p_exp = stats.pearsonr(half_lives, T1s_exp)
    pear_R_ini, p_ini = stats.pearsonr(half_lives, lambdas)
    x = np.arange(0, 1001, 1)
    y = x

    # Display scatter plot data.
    # style.use('ggplot')
    plt.figure(figsize=(9, 7))
    # plt.title('Mean Translation Time' + r' $\langle T \rangle$' +
    #          ' vs. mRNA Half-Life' + r' $t_{1/2}$', fontsize=21)
    plt.xlabel('mRNA Half-Life' + r' $t_{1/2}$' + ' (s)', fontsize=22)
    plt.ylabel('Mean Translation Time' + r' $\langle T \rangle$' + ' (s)', fontsize=22)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.scatter(gene_data["t_1/2"], gene_data["total_time"], marker='.', label="Simulation Data")
    plt.plot(half_lives, m_exp * half_lives + c_exp, 'g',
             label="Fit |r = %6.2e, p = %6.2e " % (pear_R_exp, p_exp))
    #plt.plot(x, y)
    # plt.ylim(0,300)
    # plt.grid()
    plt.legend(prop={'size': 18})

    plt.savefig("overall_gene_plot.png")
    plt.show()

    # style.use('ggplot')
    plt.figure(figsize=(8, 6))
    plt.title('Initiation Time' + r' $\lambda$' + ' vs. mRNA Half-Life' + r' $t_{1/2}$', fontsize=21)
    plt.xlabel('mRNA Half-Life' + r' $t_{1/2}$' + ' (s)', fontsize=16)
    plt.ylabel('Initiation Time' + r' $\lambda$' + ' (s)', fontsize=16)
    plt.scatter(gene_data["t_1/2"], gene_data["lambda"],
                marker='.', label="Simulation Data")
    plt.plot(half_lives, m_l * half_lives + c_l, 'tab:orange',
             label="Fit |r = %6.2e, p = %6.2e " % (pear_R_ini, p_ini))
    plt.legend(prop={'size': 16})
    plt.grid()
    plt.savefig("overall_gene_plot_ini.png")
    plt.show()


main()
