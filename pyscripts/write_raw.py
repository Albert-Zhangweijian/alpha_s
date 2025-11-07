import numpy as np


def write_raw_data(file_path, data):
    
    # write all data in uint16 format
    data = np.array(data, dtype=np.uint16)
    data.tofile(file_path)
    print(f"Raw data written to {file_path}")
    return


if __name__ == "__main__":
    # example usage
    
    file_path = "test_raw_data.bin"
    data = np.random.randint(6000, 8000, size=6400*10)  # generate some random data
    write_raw_data(file_path, data)

    