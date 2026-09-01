import matplotlib.pyplot as plt
import numpy
from os.path import join, exists
import re
import sys
plt.switch_backend('PDF')

apd90_col_idx = 13

models = ['ohara_rudy_2011', 'tentusscher_model_2006_epi', 'grandi_pasqualini_bers_2010']
#models = ['ohara_rudy_2011']

#model = 'ohara_rudy_2011'
#model = 'tentusscher_model_2006_epi'
variants = ['AZ_orig_data', 'GSK_data', 'AZ_Gintant_hERG']
variant_titles = ['Q', 'B \& Q2', 'M \& Q']

#variant = 'AZ_Gintant_hERG'
#variants = ['AZ_orig_data']
#variant = 'GSK_data'

positive_effect_cut_off = 5  # ms

drug_name_file = '../../input_data/AZ_HTS_data.txt'
file = open(drug_name_file, 'r')
drug_names = []
for line in file.readlines()[1:]:
    words = line.split()
    drug_names.append(words[0])

print drug_names

for i, name in enumerate(drug_names):
    print("{0}: {0}".format(i, name))

experimental_data_file = open('clinical_tqt_results/experimental_data.txt', 'r')
compounds = []
delta_qt_c = []
max_free_conc_micro_molar = []
experimental_data = dict()
for line in experimental_data_file.readlines()[1:]:
    drug_name, effect, conc_nM = line.split()
    experimental_data[drug_name] = [float(effect), 0.001 * float(conc_nM)]

# Select the drugs to plot here.
drugs_of_interest = drug_names


for drug_index, drug_name in enumerate(drugs_of_interest):
    print "Plotting {0}".format(drug_name)
    fig = plt.figure(drug_index, figsize=(8, 2.5))
    fig.text(0.51, 0.9, r'{0}'.format(drug_name), ha='center', va='center', fontsize=16)
    for v, variant in enumerate(variants):
        ax = plt.subplot(1, 3, v + 1)
        ax.set_title(r"{0}".format(re.sub('_', ' ', variant_titles[v])), fontsize=14)
        if v == 0:
            ax.set_ylabel(r'Change in APD$_{90}$ (ms)')
        ax.set_xlabel(r'Conc ($\mu$M)')
        ylim = [-10, 20]
        ax.set_ylim(ylim)        
        ax.set_yticklabels([r'$-10$', r'$-5$', r'$0$', r'$5$', r'$10$', r'$15$', r'$20$'])
        ax.set_xticklabels([r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$'])

        for m, model in enumerate(models):

            data_folder = join(model, variant)
            doses_for_cut_off = numpy.zeros((len(models), len(drugs_of_interest))) + 1000
            doses_95_for_cut_off = numpy.zeros((len(models), len(drugs_of_interest))) + 1000

            filename = join('../', data_folder, 'voltage_results_' + drug_name + '.dat')
            if not exists(filename):
                continue
            data = numpy.loadtxt(filename, skiprows=1)

            num_random_runs = int(max(data[:, 0]))

            for run in range(num_random_runs, -1, -1):
                row_selection = data[:, 0] == run

                doses = data[row_selection, 11]
                apd90s = data[row_selection, apd90_col_idx]
                control_apd = apd90s[0]  # Should be the same for each run, random or not.

                if run == 0:
                    # This is the original 'no noise' trace.
                    if m == 0:
                        colour_to_use = 'b'  # O'Hara
                    elif m == 1:
                        colour_to_use = 'r'  # Ten Tusscher 2006
                    elif m == 2:
                        colour_to_use = 'g'  # Grandi 2010
                    else:
                        print('Dunno what model this is!')
                        sys.exit(1)
                    ax.plot(numpy.log10(doses[1:]), apd90s[1:] - control_apd, color=colour_to_use, lw=2)

            if num_random_runs > 0:
                # Calculate two-tailed confidence interval
                concs = numpy.unique(data[:, 11])
                confidence_wanted = 95 # %

                # We want the following indices in the list of sorted random runs
                lower_idx = int(numpy.floor(((100.0 - confidence_wanted) / 2.0) * (num_random_runs / 100.0))) - 1
                upper_idx = num_random_runs - lower_idx - 2

                upper_confidence_interval = numpy.zeros([len(concs), 1])
                lower_confidence_interval = numpy.zeros([len(concs), 1])
                for conc_index in range(len(concs)):
                    this_conc_rows = data[:, 11] == concs[conc_index]
                    sorted_apd90s = sorted(data[this_conc_rows, apd90_col_idx])
                    upper_confidence_interval[conc_index] = sorted_apd90s[upper_idx] - control_apd
                    lower_confidence_interval[conc_index] = sorted_apd90s[lower_idx] - control_apd

                ax.fill_between(numpy.log10(concs[1:]),
                                upper_confidence_interval[1:, 0],
                                lower_confidence_interval[1:, 0],
                                alpha=0.2, color=colour_to_use, lw=0)

        ax.set_xlim([numpy.log10(min(doses[1:])), numpy.log10(max(doses[1:]) + 0.25 * max(doses[1:]))])
        # Plot the magic cut off line
        ax.plot(numpy.log10([doses[2], doses[-1]]), [positive_effect_cut_off, positive_effect_cut_off], 'b:', zorder=0)
        # Plot experimental data too
        ax.plot(numpy.log10([doses[2], doses[-1]]), [experimental_data[drug_name][0], experimental_data[drug_name][0]], 'k--', zorder=0)
        ax.plot(numpy.log10(experimental_data[drug_name][1]) * numpy.ones(2), ylim, 'k--', zorder=0)
        ax.plot(numpy.log10(experimental_data[drug_name][1]), experimental_data[drug_name][0], 'ro')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.subplots_adjust(top=0.75, wspace=0.25)
    plt.savefig('images/{0}.pdf'.format(drug_name), bbox_inches='tight', dpi=900, pad_inches=0.05)
