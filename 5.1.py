import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib
matplotlib.use('TkAgg')

parameters = {
    'amplitude': 1.0,
    'frequency': 2.0,
    'phase': 0.0,
    'noise_mean': 0.0,
    'noise_cov': 0.1,
    'show_noise': True,
    'filter_cutoff': 0.05
}

t = np.linspace(0, 2 * np.pi, 1000)
current_noise = np.random.normal(parameters['noise_mean'], np.sqrt(parameters['noise_cov']), size=len(t))

def harmonic(amplitude, frequency, phase):
    return amplitude * np.sin(frequency * t + phase)

def harmonic_noise(amplitude, frequency, phase, noise, show_noise):
    signal = harmonic(amplitude, frequency, phase)
    return signal + noise if show_noise else signal

def filter_signal(signal, cutoff):
    b, a = butter(N=4, Wn=cutoff, btype='low')
    return filtfilt(b, a, signal)

initial_signal = harmonic_noise(
    parameters['amplitude'],
    parameters['frequency'],
    parameters['phase'],
    current_noise,
    parameters['show_noise']
)

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.55)

signal_line, = ax.plot(t, initial_signal, label='Сигнал з шумом')
filtered_line, = ax.plot(t, filter_signal(initial_signal, parameters['filter_cutoff']),
                         label='Відфільтрований сигнал', linestyle='--')

ax.legend()

slider_axes = {
    'amp': plt.axes([0.25, 0.45, 0.65, 0.03]),
    'freq': plt.axes([0.25, 0.40, 0.65, 0.03]),
    'phase': plt.axes([0.25, 0.35, 0.65, 0.03]),
    'mean': plt.axes([0.25, 0.30, 0.65, 0.03]),
    'cov': plt.axes([0.25, 0.25, 0.65, 0.03]),
    'cutoff': plt.axes([0.25, 0.20, 0.65, 0.03])
}

sliders = {
    'amp': Slider(slider_axes['amp'], 'Амплітуда', 0.0, 5.0, valinit=parameters['amplitude']),
    'freq': Slider(slider_axes['freq'], 'Частота', 0.1, 10.0, valinit=parameters['frequency']),
    'phase': Slider(slider_axes['phase'], 'Фаза', -np.pi, np.pi, valinit=parameters['phase']),
    'mean': Slider(slider_axes['mean'], 'Середнє шуму', -1.0, 1.0, valinit=parameters['noise_mean']),
    'cov': Slider(slider_axes['cov'], 'Дисперсія шуму', 0.01, 1.0, valinit=parameters['noise_cov']),
    'cutoff': Slider(slider_axes['cutoff'], 'См. частота фільтру', 0.01, 0.5, valinit=parameters['filter_cutoff'])
}

check_ax = plt.axes([0.025, 0.7, 0.20, 0.15])
check = CheckButtons(check_ax, ['Показати шум'], [parameters['show_noise']])

reset_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_ax, 'Reset')

def update(_=None):
    amp = sliders['amp'].val
    freq = sliders['freq'].val
    phase = sliders['phase'].val
    mean = sliders['mean'].val
    cov = sliders['cov'].val
    cutoff = sliders['cutoff'].val
    show_noise = check.get_status()[0]

    global current_noise

    if mean != parameters['noise_mean'] or cov != parameters['noise_cov']:
        current_noise = np.random.normal(mean, np.sqrt(cov), size=len(t))
        parameters['noise_mean'] = mean
        parameters['noise_cov'] = cov

    parameters['filter_cutoff'] = cutoff

    signal = harmonic_noise(amp, freq, phase, current_noise, show_noise)

    signal_line.set_ydata(signal)
    filtered_line.set_ydata(filter_signal(signal, cutoff))
    fig.canvas.draw_idle()

for slider in sliders.values():
    slider.on_changed(update)

check.on_clicked(update)

def reset(event):
    for slider in sliders.values():
        slider.reset()
    if not check.get_status()[0]:
        check.set_active(0)
    update()

reset_button.on_clicked(reset)

plt.show()
