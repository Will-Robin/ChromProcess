import numpy as np


def adjacent_average(data: np.ndarray, window: int) -> np.ndarray:
    kernel_size = window
    kernel = np.ones(kernel_size) / kernel_size
    data_convolved = np.convolve(data, kernel, mode="same")
    return data_convolved
