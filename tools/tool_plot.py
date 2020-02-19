import pymzml
import numpy as np


def plot_spec_list(spec_list, mz_range, noise, output_folder):
    for spec in spec_list:
        temp_spec = spec.peaks_cut_noise(mod=noise)
        max_intensity = np.max(temp_spec[:, 1])
        temp_spec[:, 1] = temp_spec[:, 1] / max_intensity * 1000.0
        p = pymzml.plot.Factory()
        p.new_plot()
        p.add(
            temp_spec.tolist(),
            color=(0, 0, 0),
            style='sticks',
            name='peaks',
            mz_range=mz_range,
            opacity=1.0,
            title=spec.spec_id
        )
        if spec.prec_type == '[M+H]+':
            prec_type = 'H'
        elif spec.prec_type == 'M+Na]+':
            prec_type = 'Na'
        else:
            prec_type = 'NoneType'
        filename = output_folder + '{}_{}_{}_{}_{}.html'.format(spec.snapeaks_id, spec.instrument_type,
                                                                prec_type,
                                                                spec.collision_energy, spec.spec_id)
        p.save(filename=filename, mz_range=mz_range)
        print('Plotted file: {0}'.format(filename))


def plot_spec(spec, mz_range, output_file):
    temp_spec = spec.peaks_cut_noise(mod='dynamic')
    max_intensity = np.max(temp_spec[:, 1])

    temp_spec[:, 1] = temp_spec[:, 1] / max_intensity * 1000.0
    p = pymzml.plot.Factory()
    p.new_plot()
    p.add(
        temp_spec.tolist(),
        color=(0, 0, 0),
        style='sticks',
        name='peaks',
        mz_range=mz_range,
        opacity=1.0,
        title=spec.spec_id
    )

    p.save(filename=output_file, mz_range=mz_range)
    print('Plotted file: {0}'.format(output_file))
