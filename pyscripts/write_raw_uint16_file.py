import numpy as np


def write_raw_uint16_file(filepath, data):
    # write all data in uint16 format
    data = np.array(data, dtype=np.uint16)
    data.tofile(filepath)
    print(f"Raw data written to {filepath}")
    return


if __name__ == "__main__":
    # example usage
    
    file_path = "test_raw_data.bin"
    data = np.random.randint(6000, 8000, size=6400*10)  # generate some random data
    write_raw_uint16_file(file_path, data)

    