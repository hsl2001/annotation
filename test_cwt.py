#!/usr/bin/env python3
"""Validate Morlet CWT against analytic responses and library references.

Run with: uv run --python 3.12 --no-project --with scipy --with PyWavelets --with fcwt --with matplotlib python test_cwt.py
"""

import argparse
import ctypes
import csv
import importlib.util
import pathlib
import subprocess
import sys
import tempfile
import functools
import hashlib
import os
import shlex

import numpy as np
import pywt
from scipy.signal import convolve
from scipy.integrate import quad
from plot_cwt import load_matrix


SCALES = (4, 5, 6, 7, 8, 9)
CMOR = "cmor2-0.954929658551372"
CMOR_CENTER = 0.954929658551372
BASE_VALUES = {"A": 1 + 0j, "C": 1j, "G": -1j, "T": -1 + 0j}


def encode(sequence):
    return np.array([BASE_VALUES.get(base.upper(), 0j) for base in sequence], dtype=np.complex128)


def reverse_complement(sequence):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base.upper(), "N") for base in sequence[::-1])


def mother(position):
    normalization = (np.sqrt(np.pi) * (1 + np.exp(-36) - 2 * np.exp(-27))) ** -0.5
    return normalization * np.exp(-0.5 * position**2) * (np.exp(6j * position) - np.exp(-18))


def interval_integral(lower, upper):
    lower, upper = max(lower, -12), min(upper, 12)
    if lower >= upper:
        return 0j
    real = quad(lambda position: mother(position).real, lower, upper, epsabs=1e-13)[0]
    imag = quad(lambda position: -mother(position).imag, lower, upper, epsabs=1e-13)[0]
    return real + 1j * imag


@functools.lru_cache(maxsize=None)
def morlet_kernel(scale):
    radius = int(np.ceil(8 * scale + 0.5))
    return np.array([np.sqrt(scale) * interval_integral((position - 0.5) / scale,
                                                       (position + 0.5) / scale)
                     for position in range(-radius, radius + 1)])


def scipy_reference(sequence, width):
    return convolve(encode(sequence), morlet_kernel(width)[::-1], mode="same", method="direct")


def load_fcwt():
    spec = importlib.util.find_spec("fcwt")
    if spec is None or not spec.submodule_search_locations:
        raise RuntimeError("fcwt is required; install it with the command in this script's docstring")

    package = pathlib.Path(next(iter(spec.submodule_search_locations)))
    if sys.platform.startswith("linux"):
        for library in ("libfftw3fl.so", "libfftw3f_ompl.so"):
            path = package / library
            if path.is_file():
                ctypes.CDLL(str(path), mode=ctypes.RTLD_GLOBAL)

    import fcwt

    return fcwt


def verify_pywavelets_kernels():
    wavelet = pywt.ContinuousWavelet(CMOR)
    samples, sample_x = wavelet.wavefun(level=14)
    normalization = (np.sqrt(np.pi) * (1 + np.exp(-36) - 2 * np.exp(-27))) ** -0.5
    sampled = normalization * (np.sqrt(2 * np.pi) * samples - np.exp(-18 - sample_x**2 / 2))
    expected = mother(sample_x)
    np.testing.assert_allclose(sampled, expected, rtol=2e-6, atol=2e-6)
    return float(np.max(np.abs(sampled - expected)))


def verify_pywavelets_cwt(sequence):
    signal = encode(sequence)
    scales = np.asarray(SCALES, dtype=float)
    actual, _ = pywt.cwt(signal, scales, CMOR, method="conv", precision=16)
    normalization = (np.sqrt(np.pi) * (1 + np.exp(-36) - 2 * np.exp(-27))) ** -0.5
    actual *= normalization * np.sqrt(2 * np.pi)
    relative_errors = []
    correlations = []

    for index, width in enumerate(SCALES):
        radius = int(np.ceil(8 * width + 0.5))
        kernel = np.array([np.sqrt(width) * interval_integral(position / width, (position + 1) / width)
                           for position in range(-radius, radius + 1)])
        expected = convolve(signal, kernel[::-1], mode="same", method="direct")
        observed = actual[index]
        relative_errors.append(float(np.linalg.norm(observed - expected) / np.linalg.norm(expected)))
        correlations.append(float(abs(np.vdot(observed, expected)) /
                                  (np.linalg.norm(observed) * np.linalg.norm(expected))))

    if max(relative_errors) > 0.005 or min(correlations) < 0.999:
        raise AssertionError("PyWavelets integral-discretized CWT differs beyond tolerance")
    return max(relative_errors), min(correlations)


def verify_fcwt(fcwt):
    rng = np.random.default_rng(7)
    signal = (rng.normal(size=4096) + 1j * rng.normal(size=4096)).astype(np.complex64)
    sample_rate = 256
    wavelet = fcwt.Morlet(2.0)
    scale_grid = fcwt.Scales(wavelet, fcwt.FCWT_LOGSCALES, sample_rate, 8.0, 64.0, 5)
    scales = np.zeros(5, dtype=np.float32)
    frequencies = np.zeros(5, dtype=np.float32)
    scale_grid.getScales(scales)
    scale_grid.getFrequencies(frequencies)
    if not np.all(np.diff(scales) > 0) or not np.allclose(scales * frequencies, sample_rate):
        raise AssertionError("Unexpected FCWT scale ordering or frequency mapping")

    observed = np.zeros((len(scales), signal.size), dtype=np.complex64)
    fcwt.FCWT(wavelet, 1, False, True).ccwt(signal, scale_grid, observed)
    reference, _ = pywt.cwt(signal, scales, "cmor8-1.0", method="conv")

    crop = int(np.ceil(9 * scales[-1]))
    maximum_raw_error = 0.0
    maximum_adjusted_error = 0.0
    minimum_correlation = 1.0
    for index in range(len(scales)):
        actual = observed[index, crop:-crop]
        expected = reference[index, crop:-crop]
        gain = np.vdot(expected, actual) / np.vdot(expected, expected)
        raw_error = float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))
        adjusted_error = float(np.linalg.norm(actual - gain * expected) / np.linalg.norm(actual))
        correlation = float(abs(np.vdot(expected, actual)) /
                            (np.linalg.norm(expected) * np.linalg.norm(actual)))
        maximum_raw_error = max(maximum_raw_error, raw_error)
        maximum_adjusted_error = max(maximum_adjusted_error, adjusted_error)
        minimum_correlation = min(minimum_correlation, correlation)

    np.testing.assert_array_equal(np.isfinite(observed), True)
    if maximum_adjusted_error > 0.12 or minimum_correlation < 0.98:
        raise AssertionError("FCWT and PyWavelets disagree beyond the continuous-Morlet tolerance")

    return maximum_raw_error, maximum_adjusted_error, minimum_correlation


def verify_analytic():
    normalization = (np.sqrt(np.pi) * (1 + np.exp(-36) - 2 * np.exp(-27))) ** -0.5

    integral = quad(lambda position: mother(position).real, -12, 12, epsabs=1e-12)[0]
    energy = quad(lambda position: abs(mother(position))**2, -12, 12, epsabs=1e-12)[0]
    np.testing.assert_allclose(integral, 0, atol=1e-12)
    np.testing.assert_allclose(energy, 1, atol=1e-12)
    maximum_error = 0.0
    for scale in (1, 2, 3, *SCALES, 16, 64):
        np.testing.assert_allclose(morlet_kernel(scale).sum(), 0, atol=1e-12)
        for frequency in (0, 3, 6, 9, -6):
            expected = normalization * np.sqrt(2 * np.pi * scale) * (
                np.exp(-0.5 * (frequency - 6)**2) -
                np.exp(-18 - 0.5 * frequency**2))
            actual = np.sqrt(scale) * (quad(lambda position: (np.exp(1j * frequency * position) *
                np.conj(mother(position))).real, -12, 12, epsabs=1e-12)[0] +
                1j * quad(lambda position: (np.exp(1j * frequency * position) *
                np.conj(mother(position))).imag, -12, 12, epsabs=1e-12)[0])
            maximum_error = max(maximum_error, float(abs(actual - expected)))
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    return maximum_error


def verify_export(binary, scales=SCALES):
    rng = np.random.default_rng(42)
    sequences = {
        "long": "".join(rng.choice(list("ACGTNacgtX"), size=1037)),
        "short": "ACGTNacgta",
        "tiny": "A",
        "constant": "A" * 2048,
        "sinusoid": "ACTG" * 512,
        "impulse": "N" * 128 + "A" + "N" * 128,
    }

    with tempfile.TemporaryDirectory() as temporary:
        root = pathlib.Path(temporary)
        fasta = root / "test.fa"
        output = root / "cwt"
        fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in sequences.items()),
                         encoding="ascii")
        subprocess.run([str(binary), str(fasta), str(output)], check=True)

        with (output / "contigs.tsv").open() as stream:
            if stream.readline().strip() != "# anno-cwt-v2":
                raise AssertionError("Unexpected CWT index version")
            fields = stream.readline().rstrip().split("\t")
            widths = [float(value) for value in fields[1].split(",")]
            if fields[0] != "# scales" or widths != list(scales):
                raise AssertionError(f"Unexpected CWT scales: {widths}")
            header = stream.readline().rstrip().split("\t")
            if header != ["seqid", "length", "plus_offset", "minus_offset"]:
                raise AssertionError(f"Unexpected CWT index columns: {header}")
            metadata = list(csv.DictReader(stream, fieldnames=header, delimiter="\t"))

        rows = sum(2 * len(sequence) for sequence in sequences.values())
        raw = np.fromfile(output / "matrix.bin", dtype="<f4")
        if raw.size != rows * len(scales) * 2:
            raise AssertionError(f"Unexpected matrix size: {raw.size} float32 values")
        observed = raw.reshape(rows, len(scales), 2)
        observed = observed[..., 0] + 1j * observed[..., 1]
        loaded, contigs, loaded_scales = load_matrix(output)
        if loaded_scales != list(scales) or any(item["parameter"] != "Scale" for item in contigs.values()):
            raise AssertionError("Plot loader misinterprets CWT scales")
        np.testing.assert_array_equal(loaded, observed)

        maximum_error = 0.0
        checked_coefficients = 0
        expected_offset = 0
        for item in metadata:
            name = item["seqid"]
            sequence = sequences[name]
            length = int(item["length"])
            plus = int(item["plus_offset"])
            minus = int(item["minus_offset"])
            if length != len(sequence) or plus != expected_offset or minus != plus + length:
                raise AssertionError(f"Invalid offsets for {name}")

            for offset, strand_sequence in ((plus, sequence), (minus, reverse_complement(sequence))):
                for scale, width in enumerate(scales):
                    expected = scipy_reference(strand_sequence, width)
                    actual = observed[offset:offset + length, scale]
                    error = float(np.max(np.abs(actual - expected)))
                    maximum_error = max(maximum_error, error)
                    checked_coefficients += length
                    np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=2e-7)
                    if name == "constant":
                        radius = int(np.ceil(8 * width + 0.5))
                        np.testing.assert_allclose(actual[radius:-radius], 0, atol=2e-7)
                    elif name == "sinusoid":
                        radius = int(np.ceil(8 * width + 0.5))
                        positions = np.arange(radius, length - radius)
                        frequency = np.pi / 2
                        if offset == minus:
                            frequency = -frequency
                        normalization = (np.sqrt(np.pi) * (1 + np.exp(-36) - 2 * np.exp(-27))) ** -0.5
                        harmonic_count = int(np.ceil(18 / (2 * np.pi * width))) + 2
                        frequencies = frequency + 2 * np.pi * np.arange(-harmonic_count, harmonic_count + 1)
                        response = normalization * np.sqrt(2 * np.pi * width) * np.sum(
                            np.sinc(frequencies / (2 * np.pi)) * (
                            np.exp(-0.5 * (width * frequencies - 6)**2) -
                            np.exp(-18 - 0.5 * (width * frequencies)**2)))
                        initial = encode(strand_sequence)[0]
                        analytic = initial * np.exp(1j * frequency * positions) * response
                        np.testing.assert_allclose(actual[positions], analytic, rtol=2e-6, atol=2e-7)
                    elif name == "impulse":
                        center_value = encode(strand_sequence)[128]
                        positions = np.arange(length)
                        analytic = center_value * np.sqrt(width) * np.array([
                            interval_integral((127.5 - position) / width, (128.5 - position) / width)
                            for position in positions])
                        np.testing.assert_allclose(actual, analytic, rtol=2e-6, atol=2e-7)
            expected_offset += 2 * length

        if expected_offset != rows:
            raise AssertionError("CWT index does not cover the full matrix")

    return checked_coefficients, maximum_error, hashlib.sha256(raw.tobytes()).hexdigest()


def verify_configurations():
    source = pathlib.Path(__file__).resolve().with_name("anno_cwt.c")
    scales = (1, 2, 3, 7, 16, 64)
    digests = []
    total = 0
    maximum_error = 0.0
    with tempfile.TemporaryDirectory() as temporary:
        binary = pathlib.Path(temporary) / "anno_cwt"
        compiler = shlex.split(os.environ.get("CC", "cc"))
        for capacity in (0, 1, 8, 145, 4096):
            subprocess.run([*compiler, "-O2", "-std=c11", "-Wall", "-Wextra",
                            "-DCWT_SCALES=" + ",".join(map(str, scales)),
                            f"-DMAX_WAVE_SIZE={capacity}", str(source), "-lz", "-lm", "-o", str(binary)], check=True)
            checked, error, digest = verify_export(binary, scales)
            total += checked
            maximum_error = max(maximum_error, error)
            digests.append(digest)
        if len(set(digests)) != 1:
            raise AssertionError("MAX_WAVE_SIZE changes mathematical results")
        duplicate_scales = (32, 3, 3, 1)
        subprocess.run([*compiler, "-O2", "-std=c11", "-DCWT_SCALES=" + ",".join(map(str, duplicate_scales)),
                        "-DMAX_WAVE_SIZE=2", str(source), "-lz", "-lm", "-o", str(binary)], check=True)
        checked, error, _ = verify_export(binary, duplicate_scales)
        total += checked
        maximum_error = max(maximum_error, error)
        for invalid in ("0", "-1", "0.5", "NAN", "INFINITY", "1e300"):
            subprocess.run([*compiler, "-O2", "-std=c11", f"-DCWT_SCALES={invalid}", str(source),
                            "-lz", "-lm", "-o", str(binary)], check=True)
            result = subprocess.run([str(binary), "unused.fa", str(binary.parent / "out")], capture_output=True, text=True)
            if result.returncode == 0 or "CWT scale" not in result.stderr:
                raise AssertionError(f"Invalid scale not rejected: {invalid}")
    return total, maximum_error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=pathlib.Path,
                        default=pathlib.Path(__file__).resolve().with_name("anno_cwt"))
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"CWT executable not found: {binary}; run make anno_cwt first")

    analytic_error = verify_analytic()
    kernel_delta = verify_pywavelets_kernels()
    checked, coefficient_error, _ = verify_export(binary)
    configured, configuration_error = verify_configurations()
    sequence = "".join(np.random.default_rng(23).choice(list("ACGTN"), size=257))
    pywt_error, pywt_correlation = verify_pywavelets_cwt(sequence)
    fcwt = load_fcwt()
    fcwt_raw_error, fcwt_adjusted_error, fcwt_correlation = verify_fcwt(fcwt)
    print(f"PASS analytic Morlet: zero mean, unit energy, Fourier response; max abs error={analytic_error:.3g}")
    print(f"PASS PyWavelets cmor kernel samples: max abs error={kernel_delta:.3g}")
    print(f"PASS SciPy convolution: {checked:,} complex coefficients; max abs error={coefficient_error:.3g}")
    print(f"PASS arbitrary scales/capacities: {configured:,} coefficients; max abs error={configuration_error:.3g}; identical capacity-independent results; invalid scales rejected")
    print(f"PASS PyWavelets CWT (known normalization, b=n-0.5 alignment, no fitted gain): max relative error={pywt_error:.3g}; min correlation={pywt_correlation:.3g}")
    print(f"APPROX FCWT vs PyWavelets cmor8-1.0: max raw relative error={fcwt_raw_error:.3g}; "
          f"gain-adjusted residual={fcwt_adjusted_error:.3g}; min correlation={fcwt_correlation:.3g}")
    print("Covered both strands, DC/sinusoid/box-pulse continuous responses, edge padding, ambiguity bases, all scales, and the 1024-base chunk boundary.")


if __name__ == "__main__":
    main()