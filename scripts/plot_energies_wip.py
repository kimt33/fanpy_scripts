import csv
from matplotlib import pyplot as plt
import numpy as np
import sys
sys.path.append('../')
from basic_units import radians
import itertools as it


markers = np.array(['o', '<', '>', 'v', '1', '2', '3', '4', '+', 'x', '.'])
colors = np.array(['tab:pink', 'tab:brown', 'tab:gray', 'tab:red', 'tab:blue', 'tab:orange', 'tab:green',
          'tab:purple', 'tab:olive', 'tab:cyan', 'black'])
sizes = np.array([5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 5])
plot_label = {
}

system = "h8"
units = "Angstrom"

basis_sets = ['sto-6g']

# NOTE: untab and remove
plot_labels = {
    'H8_stretch_apg_100steps_lstsq' : {
        "fci" : "FCI",
        "hf" : "HF",
        "apg" : "APG",
        "apg_fanpt_order1_100steps_lstsq" : "APG_FanPT_1",
        "apg_fanpt_order2_100steps_lstsq" : "APG_FanPT_2",
        "apg_fanpt_order3_100steps_lstsq" : "APG_FanPT_3",
        "apg_fanpt_order4_100steps_lstsq" : "APG_FanPT_4",
    },
    'H8_stretch_apg_10steps_lstsq' : {
        "fci" : "FCI",
        "hf" : "HF",
        "apg" : "APG",
        "apg_fanpt_order1_10steps_lstsq" : "APG_FanPT_1",
        "apg_fanpt_order2_10steps_lstsq" : "APG_FanPT_2",
        "apg_fanpt_order3_10steps_lstsq" : "APG_FanPT_3",
        "apg_fanpt_order4_10steps_lstsq" : "APG_FanPT_4",
    },
    'H8_stretch_apg_100steps_pure' : {
        "fci" : "FCI",
        "hf" : "HF",
        "apg" : "APG",
        "apg_fanpt_order1_100steps_pure" : "APG_FanPT_1",
        "apg_fanpt_order2_100steps_pure" : "APG_FanPT_2",
        "apg_fanpt_order3_100steps_pure" : "APG_FanPT_3",
        "apg_fanpt_order4_100steps_pure" : "APG_FanPT_4",
    },
    'H8_stretch_apg_10steps_pure' : {
        "fci" : "FCI",
        "hf" : "HF",
        "apg" : "APG",
        "apg_fanpt_order1_10steps_pure" : "APG_FanPT_1",
        "apg_fanpt_order2_10steps_pure" : "APG_FanPT_2",
        "apg_fanpt_order3_10steps_pure" : "APG_FanPT_3",
        "apg_fanpt_order4_10steps_pure" : "APG_FanPT_4",
    },
}

# plot parameter
dpi = 400
for savename, plot_label in plot_labels.items():
    calc_types = list(plot_label.keys())

    output_table = []
    data = {}
    for basis in basis_sets:
        for calc_type in calc_types:
            filename = f'{system}_{basis}_{calc_type}_results.csv'
            input_table = []
            with open(filename, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                    input_table.append(row)

                if not output_table:
                    output_table.append([input_table[0][1], input_table[0][3].replace('Electronic', 'Total')])
                else:
                    output_table[0].append(input_table[0][3].replace('Electronic', 'Total'))
                data[calc_type] = {}

                for row in input_table[1:]:
                    if not row[3]:
                        continue
                    try:
                        data[calc_type][float(row[1])].append(float(row[2]) + float(row[3]))
                    except KeyError:
                        data[calc_type][float(row[1])] = [float(row[2]) + float(row[3])]

        for alpha, energies in data['fci'].items():
            data['fci'][alpha] = min(energies)

        print(system)
        print(data.keys())
        alphas = set()
        for calc_type, calc_type_data in data.items():
            if calc_type == 'fci':
                continue
            for alpha, energies in calc_type_data.items():
                alphas.add(alpha)
                data[calc_type][alpha] = energies[np.argmin(np.abs(np.array(energies) - np.array(data['fci'][alpha])))]

        # make graphs
        graph_title = "H8 Symmetric Stretch"

        plt.figure()

        markers_graph = iter(['o', 'o', 'o', 'o', '<', '>', '1', '2', '3', '.'])
        markers_graph = it.cycle(['o', 'v', '^', '<', '>', 's', 'p', 'P', '*'])
        colors_graph = iter(['tab:pink', 'tab:purple', 'tab:green', 'tab:orange', 'tab:blue', 'tab:red', 'black'])
        colors_graph = it.cycle(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan'])
        sizes_graph = iter([5, 5, 5, 5, 5, 5, 8, 8, 8, 5])
        sizes_graph = it.cycle([5, 5, 5, 5, 5, 5, 5, 5, 5, 5])
        for calc_type in calc_types:
            calc_type_data = data[calc_type]
            alphas = np.array(list(calc_type_data.keys()))
            energies = np.array(list(calc_type_data.values()))
            alpha_sort = np.argsort(alphas)
            xunits = None
            if calc_type == 'fci':
                plt.plot(alphas[alpha_sort], energies[alpha_sort], label=plot_label[calc_type], xunits=xunits,
                         marker='.', color='black', markersize=5)
            else:
                plt.plot(alphas[alpha_sort], energies[alpha_sort], label=plot_label[calc_type], xunits=xunits,
                         marker=next(markers_graph), color=next(colors_graph), markersize=next(sizes_graph))
        plt.title(graph_title)
        plt.legend()
        plt.xlabel(f"Bond Length (Angstrom)")
        plt.ylabel("Energy (Hartree)")
        plt.tight_layout()
        plt.savefig(f'./{savename}.png', transparent=True, dpi=dpi)
        
        plt.figure()
        fci_alphas = np.array(list(data['fci'].keys()))
        fci_energies = np.array(list(data['fci'].values()))
        fci_energies = fci_energies[np.argsort(fci_alphas)]
        markers_graph = iter(['o', 'o', 'o', 'o', '<', '>', '1', '2', '3', '.'])
        markers_graph = it.cycle(['o', 'v', '^', '<', '>', 's', 'p', 'P', '*'])
        colors_graph = iter(['tab:pink', 'tab:purple', 'tab:green', 'tab:orange', 'tab:blue', 'tab:red', 'black'])
        colors_graph = it.cycle(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan'])
        sizes_graph = iter([5, 5, 5, 5, 5, 5, 8, 8, 8, 5])
        sizes_graph = it.cycle([5, 5, 5, 5, 5, 5, 5, 5, 5, 5])
        for calc_type in calc_types:
            calc_type_data = data[calc_type]
            alphas = np.array(list(calc_type_data.keys()))
            energies = np.array(list(calc_type_data.values()))
            alpha_sort = np.argsort(alphas)
            xunits = None
            if calc_type == 'fci':
                plt.plot(alphas[alpha_sort], energies[alpha_sort] - fci_energies, label=plot_label[calc_type], xunits=xunits,
                         marker='.', color='black', markersize=5)
            else:
                plt.plot(alphas[alpha_sort], energies[alpha_sort] - fci_energies, label=plot_label[calc_type], xunits=xunits,
                         marker=next(markers_graph), color=next(colors_graph), markersize=next(sizes_graph))
        plt.title(graph_title + " FCI Difference")
        plt.legend()
        plt.xlabel(f"Bond Length (Angstrom)")
        plt.ylabel("Energy Difference (Hartree)")

        plt.tight_layout()
        # XXXX
        plt.savefig(f'./{savename}_fci_diff.png', transparent=True, dpi=dpi)
        #plt.savefig('./graphs/' + system.replace('.', '') + '_fci_diff.png', transparent=True)
        #plt.savefig('./graphs/' + system.replace('.', '') + '_senzero_fci_diff.png', transparent=True)
        #plt.show()
