import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import csv

def get_correlation(dict, x_label, y_label):
    x_data = np.array(dict["x_label"])
    y_data = np.array(dict["y_label"])
    pear_r, p_val = stats.pearsonr(x_data, y_data)
    return pear_r, p_val

def get_fit(dict, x_label, y_label):
    x_data = np.array(dict["x_label"])
    y_data = np.array(dict["y_label"])
    M = np.vstack([x_data, np.ones(len(x_data))]).T
    m, c = np.linalg.lstsq(M, y_data)[0]
    return m, c


def plotter(dict, gradient, intercept, pear_r, p_val, x_label, y_label, xaxis_label, yaxis_label, title, file_name):
    # Display scatter plot data.
    plt.figure(figsize=(10, 8))
    plt.title(str(title), fontsize=20)
    plt.xlabel(str(xaxis_label), fontsize=15)
    plt.ylabel(str(yaxis_label), fontsize=15)
    plt.scatter(dict[str(x_label)], dict[str(y_label)],
                marker='o')
    x_data = np.array(dict[str(x_label)])
    plt.plot(x_data, gradient*x_data + intercept, 'r',
            label="Fit |%6s, r = %6.2e" % ('red', pear_r, p_val))
    plt.legend()
    plt.savefig(str(file_name))
    plt.show()

def main():
    # Open file.
    with open("unstable_genes.csv", "r") as f:
        reader = csv.reader(f)
        # Skip first line.
        next(f)
        # Dictionary for plotting.
        gene_dict = {}
        check = False
        # Assigning values from CSV.
        for row in reader:
            key = row[0]
            # Duplicate row handling.
            if key in gene_dict:
                pass
            # Making data floats.
            if check == False:
                gene_dict[key] = row[1:]
                check = True
            else:
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

main()
