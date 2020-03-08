import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
import scipy.stats as stats
import csv

def get_fit(dict, x_label, y_label):
    x_data = np.array(dict[x_label], dtype=np.float64)
    y_data = np.array(dict[y_label], dtype=np.float64)
    M = np.vstack([x_data, np.ones(len(x_data))]).T
    m, c = np.linalg.lstsq(M, y_data, rcond=None)[0]
    return m, c

def get_correlation(dict, x_label, y_label):
    x_data = np.array(dict[x_label], dtype=np.float64)
    y_data = np.array(dict[y_label], dtype=np.float64)
    pear_r, p_val = stats.pearsonr(x_data, y_data)
    return pear_r, p_val

def plotter(dict, gradient, intercept, pear_r, p_val, x_label, y_label, xaxis_label, yaxis_label, title, file_name):
    # Display scatter plot data.
    style.use('seaborn-poster')  
    style.use('ggplot')
    plt.title(title, fontsize=20)
    plt.xlabel(xaxis_label, fontsize=15)
    plt.ylabel(yaxis_label, fontsize=15)
    plt.scatter(dict[x_label], dict[y_label], marker='.', label = "Data")
    x_data = np.array(dict[x_label], dtype=np.float64)
    plt.plot(x_data, gradient*x_data + intercept, 'b',
             label="Fit | r = %6.2e, p = %6.2e" % (pear_r, p_val))
    plt.legend()
    plt.savefig('unstable_plots/' + file_name)
    plt.show()

def main():
    # Open file.
    with open("unstable_genes.csv", "r") as f:
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
    gene_data = {
                "T1_exp": [], "T1_ana": [], "t_1/2": [], "lambda": [],
                "CDS": [], "alpha": [], "label": []
                }
    for label, coord in gene_dict.items():
        gene_data["T1_exp"].append(coord[0])
        gene_data["T1_ana"].append(coord[1])
        gene_data["t_1/2"].append(coord[2])
        gene_data["lambda"].append(coord[3])
        gene_data["CDS"].append(coord[4])
        gene_data["alpha"].append(coord[5])
        gene_data["label"].append(label)

    # Plotting.
    # Translation Time vs. mRNA Half-Life.
    m_1, c_1 = get_fit(gene_data, "t_1/2", "T1_exp")
    r_1, p_1 = get_correlation(gene_data, "t_1/2", "T1_exp")
    plotter(gene_data, m_1, c_1, r_1, p_1, "t_1/2",
            "T1_exp", 'mRNA Half-Life (s)', 'Mean Translation Time (s)', 'Translation Time vs. mRNA Half-Life',
            "unstable_T1_half.png")

    # Initiation Time vs. Half-life.
    m_2, c_2 = get_fit(gene_data, "t_1/2", "lambda")
    r_2, p_2 = get_correlation(gene_data, "t_1/2", "lambda")
    plotter(gene_data, m_2, c_2, r_2, p_2, "t_1/2",
            "lambda", 'mRNA Half-Life (s)', 'Initiation Time ' +
            r'[$\lambda$] ' + '(s)', 'Initiation Time vs. mRNA Half-Life',
            "unstable_lambda_half.png")

    # CDS vs. Half-life
    m_3, c_3 = get_fit(gene_data, "t_1/2", "CDS")
    r_3, p_3 = get_correlation(gene_data, "t_1/2", "CDS")
    plotter(gene_data, m_3, c_3, r_3, p_3, "t_1/2",
            "CDS", 'mRNA Half-Life (s)', 'CDS Length (codons)', 'CDS vs. mRNA Half-Life',
            "unstable_cds_half.png")
    
    # Initiation Rate vs. Half-life.
    m_4, c_4 = get_fit(gene_data, "t_1/2", "alpha")
    r_4, p_4 = get_correlation(gene_data, "t_1/2", "alpha")
    plotter(gene_data, m_4, c_4, r_4, p_4, "t_1/2",
            "alpha", 'mRNA Half-Life (s)', r"[$\alpha$] $s^{-1}$", 'Initiation Rate vs. mRNA Half-Life',
            "unstable_alpha_half.png")

    # Translation Time vs. CDS
    m_5, c_5 = get_fit(gene_data, "CDS", "T1_exp")
    r_5, p_5 = get_correlation(gene_data, "CDS", "T1_exp")
    plotter(gene_data, m_5, c_5, r_5, p_5, "CDS",
            "T1_exp", 'CDS Length (codons)', 'Mean Translation Time (s)', 'Translation Time vs. CDS',
            "unstable_t1_cds.png")

    # Translation Time vs. Initiation Time.
    m_6, c_6 = get_fit(gene_data, "lambda", "T1_exp")
    r_6, p_6 = get_correlation(gene_data, "lambda", "T1_exp")
    plotter(gene_data, m_6, c_6, r_6, p_6, "lambda",
            "T1_exp", 'Initiation Time ' +
            r'[$\lambda$] ', 'Mean Translation Time (s)', 'Translation Time vs. Initiation Time',
            "unstable_t1_lambda.png")
        
    # Translation Time vs. Initiation Rate.
    m_7, c_7 = get_fit(gene_data, "alpha", "T1_exp")
    r_7, p_7 = get_correlation(gene_data, "alpha", "T1_exp")
    plotter(gene_data, m_7, c_7, r_7, p_7, "alpha",
            "T1_exp", 'Initiation Rate ' +
            r'[$\alpha$] $s^{-1}$', 'Mean Translation Time (s)', 'Translation Time vs. Initiation Rate',
            "unstable_t1_alpha.png")
    
    # Translation Time vs. CDS
    m_8, c_8 = get_fit(gene_data, "CDS", "alpha")
    r_8, p_8 = get_correlation(gene_data, "CDS", "alpha")
    plotter(gene_data, m_8, c_8, r_8, p_8, "CDS",
            "alpha", 'CDS Length (codons)', r'[$\alpha$] $s^{-1}$', 'Initiation Rate vs. CDS',
            "unstable_alpha_cds.png")
    
main()
