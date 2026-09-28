"""OFDM demo: modulate, pass through a simple channel, then acquire and track in the receiver.

    bazel run //:ofdm_demo                 # interactive windows (needs a display)
    bazel run //:ofdm_demo -- --png /tmp   # write the figures as PNGs instead

Transmitter: 512-subcarrier OFDM with a 1/8 cyclic prefix at 10 MHz. Every 40 symbols, three
BPSK Gold-code reference symbols (120 centre subcarriers), then QPSK and 16-QAM data symbols.
Channel: unknown start delay, noise, and Doppler (carrier offset + time compression).
Receiver: coarse correlation over everything, then symbol-by-symbol demodulation, then a
timing tracking loop that re-finds the correlation peak at each reference burst.

Things to try in the grid window: put the pointer on a QPSK or 16-QAM row, press x for that
symbol's cut, then 5 for its constellation (Esc to go back).
"""

import argparse
import os

import numpy as np

import signal_sniper_plot_py as ssp

# ── Parameters ──────────────────────────────────────────────────────────────────
FS = 10e6              # sample rate
N = 512                # FFT size / subcarriers
CP = N // 8            # cyclic prefix
SYM = N + CP           # samples per symbol
NUM_SYMBOLS = 100
PERIOD = 40            # reference burst every 40 symbols
NUM_REF = 3            # reference symbols per burst (one per Gold code)
REF_BINS = 120         # subcarriers used by each reference
DATA_BINS = 240        # subcarriers used by data symbols
TRACK_BINS = 128       # centre subcarriers used for the tracking correlation

SNR_DB = 15.0
CFO_HZ = 2_000.0       # Doppler carrier offset (~0.1 subcarrier)
DRIFT_PPM = 200.0      # Doppler time compression; exaggerated so it shows within 100 symbols
LEAD_SAMPLES = 1_234   # unknown delay before the first symbol

rng = np.random.default_rng(7)

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--png", metavar="DIR", help="write the figures as PNGs to DIR instead of opening windows")
args = parser.parse_args()


def centre(width):
    """Bin indices of `width` subcarriers centred on DC (bins in fftshifted order)."""
    return np.arange(N // 2 - width // 2, N // 2 + width // 2)


# ── Transmitter ─────────────────────────────────────────────────────────────────
def m_sequence(taps, degree=7):
    """One period of a maximal-length sequence from a Fibonacci LFSR."""
    state = [1] * degree
    out = []
    for _ in range(2 ** degree - 1):
        out.append(state[-1])
        fb = 0
        for t in taps:
            fb ^= state[t - 1]
        state = [fb] + state[:-1]
    return np.array(out, dtype=np.int8)


def gold_codes(count, length):
    """`count` Gold codes from the degree-7 preferred pair x^7+x^3+1 and x^7+x^3+x^2+x+1."""
    a, b = m_sequence([7, 3]), m_sequence([7, 3, 2, 1])
    return [(a ^ np.roll(b, k))[:length] for k in range(count)]


def qam(bits_per_symbol, count):
    """Random unit-power QPSK (2) or 16-QAM (4) symbols."""
    levels = 2 ** (bits_per_symbol // 2)
    axis = 2 * np.arange(levels) - (levels - 1)
    i, q = rng.choice(axis, count), rng.choice(axis, count)
    return (i + 1j * q) / np.sqrt(2 * np.mean(axis ** 2))


def ofdm_symbol(bins, values):
    """Time-domain symbol with cyclic prefix; `bins` are fftshifted subcarrier indices."""
    grid = np.zeros(N, complex)
    grid[bins] = values
    x = np.fft.ifft(np.fft.ifftshift(grid)) * np.sqrt(N)
    return np.concatenate([x[-CP:], x])


golds = gold_codes(NUM_REF, REF_BINS)
refs = [1.0 - 2.0 * g for g in golds]  # BPSK
tx_grid = np.zeros((NUM_SYMBOLS, N), complex)
for s in range(NUM_SYMBOLS):
    k = s % PERIOD
    if k < NUM_REF:
        tx_grid[s, centre(REF_BINS)] = refs[k]
    elif k < 20:
        tx_grid[s, centre(DATA_BINS)] = qam(2, DATA_BINS)
    else:
        tx_grid[s, centre(DATA_BINS)] = qam(4, DATA_BINS)
tx = np.concatenate([ofdm_symbol(np.arange(N), row) for row in tx_grid])

# ── Channel: delay, Doppler (time compression + carrier offset), noise ─────────
signal = np.concatenate([np.zeros(LEAD_SAMPLES, complex), tx, np.zeros(2 * SYM, complex)])
n = np.arange(len(signal))
t = n * (1 + DRIFT_PPM * 1e-6)  # the receiver samples a slightly compressed waveform
rx = np.interp(t, n, signal.real) + 1j * np.interp(t, n, signal.imag)
rx *= np.exp(2j * np.pi * CFO_HZ * n / FS)
noise_power = np.mean(np.abs(tx) ** 2) / 10 ** (SNR_DB / 10)
rx += np.sqrt(noise_power / 2) * (rng.standard_normal(len(rx)) + 1j * rng.standard_normal(len(rx)))

# ── Receiver 1: coarse correlation over all of the received data ────────────────
template = ofdm_symbol(centre(REF_BINS), refs[0])  # first reference symbol, with CP
corr = np.correlate(rx, template, mode="valid") / np.linalg.norm(template)  # complex
mag = np.abs(corr)                                        # peaks are found by magnitude
first = int(np.argmax(mag > 0.5 * mag.max()))             # first strong peak...
start = first + int(np.argmax(mag[first:first + SYM]))    # ...and its maximum
print(f"coarse detection: symbol 0 starts at sample {start} (true {LEAD_SAMPLES})")

# ── Receiver 2: cut into symbols, strip the CP, demodulate ──────────────────────
def demod(sym_start):
    """Subcarriers (fftshifted) of the symbol whose CP begins at `sym_start`."""
    body = rx[sym_start + CP:sym_start + SYM]
    return np.fft.fftshift(np.fft.fft(body)) / np.sqrt(N)


rx_grid = np.array([demod(start + s * SYM) for s in range(NUM_SYMBOLS)])

# ── Receiver 3: tracking loop on each reference burst ───────────────────────────
track = centre(TRACK_BINS)
ref_in_track = np.zeros((NUM_REF, TRACK_BINS))  # references on the 128-bin window (0 outside)
for j in range(NUM_REF):
    ref_in_track[j, (TRACK_BINS - REF_BINS) // 2:(TRACK_BINS + REF_BINS) // 2] = refs[j]

def true_start(s):
    """Where symbol s really starts in rx (the channel compresses time by 1 + DRIFT_PPM)."""
    return (LEAD_SAMPLES + s * SYM) / (1 + DRIFT_PPM * 1e-6)


timing = float(start)  # loop state: where symbol 0 starts, refined at each burst
rows, labels = [], []
print("burst  symbol  peak(bins)  measured error  true error  new timing   (errors in samples, + = late)")
for burst_symbol in range(0, NUM_SYMBOLS - NUM_REF + 1, PERIOD):
    errors, peaks = [], []
    for j in range(NUM_REF):
        s = burst_symbol + j
        y = demod(int(round(timing)) + s * SYM)[track]
        c = np.fft.fftshift(np.fft.ifft(y * np.conj(ref_in_track[j])))  # delay profile
        rows.append(c)
        labels.append(s)
        mag = np.abs(c)
        p = int(np.argmax(mag))
        # Parabolic interpolation around the peak; zero delay is bin TRACK_BINS // 2.
        a, b, d = mag[(p - 1) % TRACK_BINS], mag[p], mag[(p + 1) % TRACK_BINS]
        peak = p + 0.5 * (a - d) / (a - 2 * b + d) - TRACK_BINS // 2
        peaks.append(peak)
        # A window late by e samples multiplies subcarrier k by exp(+2j*pi*k*e/N), which the
        # IDFT turns into a peak at -e * TRACK_BINS / N bins.
        errors.append(-peak * N / TRACK_BINS)
    err = float(np.mean(errors))
    true_err = timing + burst_symbol * SYM - true_start(burst_symbol)
    print(f"{burst_symbol // PERIOD:5d}  {burst_symbol:6d}  {np.mean(peaks):10.2f}  "
          f"{err:15.2f}  {true_err:14.2f}  {timing - err:11.2f}")
    timing -= err  # first-order loop, gain 1
tracking = np.array(rows)

# ── Show everything ─────────────────────────────────────────────────────────────
corr_trace = ssp.Signal(corr, xdelta=1 / FS, name="xcorr with reference 0")
corr_opts = dict(cmode="mag", title="coarse correlation over the whole capture (x in seconds)")
grid_opts = dict(xstart=-N // 2, cmode="20log", zrange=(-30, 10),
                 title="received OFDM grid: subcarrier (x) vs symbol (y)")
track_opts = dict(xstart=-TRACK_BINS // 2, cmode="mag", cmap="hot",
                  title="tracking: IDFT(Y * conj(gold)), rows = 3 references per burst")
# The same delay profiles overlaid: one trace per reference symbol.
overlay = [ssp.Signal(c, xstart=-TRACK_BINS // 2,  # complex: try 2 (phase), 3/4, 5 (IR)
                      name=f"burst {s // PERIOD} gold {s % PERIOD} (sym {s})")
           for c, s in zip(rows, labels)]
overlay_opts = dict(cmode="mag",
                    title="tracking correlations overlaid (x in 128-point IDFT bins, 4 samples each)")

if args.png:
    ssp.save_png(os.path.join(args.png, "coarse_correlation.png"), corr_trace, **corr_opts)
    ssp.save_raster_png(os.path.join(args.png, "ofdm_grid.png"), rx_grid, **grid_opts)
    ssp.save_raster_png(os.path.join(args.png, "tracking.png"), tracking, **track_opts)
    ssp.save_png(os.path.join(args.png, "tracking_overlay.png"), overlay, **overlay_opts)
    print("wrote coarse_correlation.png, ofdm_grid.png, tracking.png, tracking_overlay.png to", args.png)
else:
    ssp.plot(corr_trace, **corr_opts)
    ssp.raster(rx_grid, **grid_opts)
    ssp.raster(tracking, **track_opts)
    ssp.plot(overlay, **overlay_opts)
