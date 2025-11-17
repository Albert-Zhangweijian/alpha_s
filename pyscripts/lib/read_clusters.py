import os
import struct
import time
import numpy as np

class Cluster:
    def __init__(self, frame_id, pixel_ids, adus, energies):
        self.frame_id = frame_id
        self.pixel_ids = pixel_ids
        self.adus = adus
        self.energies = energies

    @property
    def total_energy(self):
        return sum(self.energies)

    def __repr__(self):
        pixel_ids_str = "|".join([f"{pixel_id:>6d}" for pixel_id in self.pixel_ids])
        adus_str = "|".join([f"{adu:>6d}" for adu in self.adus])
        energies_str = "|".join([f"{energy:>6.2f}" for energy in self.energies])
        return f"Cluster(frame_id={self.frame_id:<d}, total energy={sum(self.energies):6.2f})keV\n\tPixel ID={pixel_ids_str}\n\tPixelADU={adus_str}\n\tPixelErg={energies_str}\n"


def read_clusters(filename, n_max=None, n_pixels=None):
    n_clusters = 0
    clusters = []
    with open(filename, "rb") as file:
        while True:
            header = file.read(8)
            if not header:
                break  

            frame_id, num_pixels = struct.unpack("ii", header)
            if n_pixels is None or num_pixels == n_pixels:
                pixel_ids = list(struct.unpack(f"{num_pixels}H", file.read(num_pixels * 2)))
                adus = list(struct.unpack(f"{num_pixels}H", file.read(num_pixels * 2)))
                energies = list(struct.unpack(f"{num_pixels}f", file.read(num_pixels * 4)))
                clusters.append(Cluster(frame_id, pixel_ids, adus, energies))
                n_clusters += 1
                if n_clusters % 1000 == 0:
                    print(f"Read {n_clusters} clusters")
            else:
                file.read(num_pixels * 8)

            if n_max is not None and n_clusters >= n_max:
                break
    print(f"Read {n_clusters} clusters")

    return clusters

def read_clusters_fast(filename):
    """
    Ultra-fast reader for fixed-size clusters (Python version of MATLAB read_clusters_parallel)
    
    File structure (repeated for each cluster):
        int32 frame_id
        int32 num_pixels
        uint16 pixel_ids[num_pixels]
        uint16 adus[num_pixels]
        float32 energies[num_pixels]
    """

    t0 = time.time()

    # --- Step 1. Read header of first cluster to determine num_pixels ---
    with open(filename, 'rb') as f:
        header = np.fromfile(f, dtype=np.int32, count=2)
    if len(header) < 2:
        raise ValueError("Corrupted file or missing header.")

    num_pixels = int(header[1])
    if num_pixels > 10:
        raise ValueError(f"This fast version assumes num_pixels ≤ 10 (current {num_pixels}).")

    # --- Step 2. Compute cluster size and total number of clusters ---
    cluster_bytes = 8 + num_pixels * (2 + 2 + 4)
    file_size = os.path.getsize(filename)
    n_clusters = file_size // cluster_bytes
    print(f"Detected {n_clusters} clusters, each {cluster_bytes} bytes ({file_size/1e6:.2f} MB total)")

    # --- Step 3. Read the whole file at once ---
    with open(filename, 'rb') as f:
        raw = np.fromfile(f, dtype=np.uint8, count=file_size)

    # --- Step 4. Parse each field ---

    # 1️⃣ frame_ids
    frame_id_indices = np.arange(0, file_size, cluster_bytes)
    frame_id_bytes = np.stack([raw[frame_id_indices + j] for j in range(4)], axis=0)
    frame_ids = np.frombuffer(frame_id_bytes.T.tobytes(), dtype=np.int32)

    # 2️⃣ pixel_ids
    pixel_start = 8  # skip frame_id(4) + num_pixels(4)
    pixel_indices = np.arange(pixel_start, file_size, cluster_bytes)
    pixel_bytes = np.stack([raw[pixel_indices + j] for j in range(2 * num_pixels)], axis=0)
    pixel_ids = np.frombuffer(pixel_bytes.T.tobytes(), dtype=np.uint16)
    pixel_ids = pixel_ids.reshape((n_clusters, num_pixels)).T.astype(np.int32)  # convert to 1-based index

    # 3️⃣ adus
    adu_start = 8 + 2 * num_pixels
    adu_indices = np.arange(adu_start, file_size, cluster_bytes)
    adu_bytes = np.stack([raw[adu_indices + j] for j in range(2 * num_pixels)], axis=0)
    adus = np.frombuffer(adu_bytes.T.tobytes(), dtype=np.uint16)
    adus = adus.reshape((n_clusters, num_pixels)).T.astype(np.int32)

    # 4️⃣ energies
    energy_start = 8 + 4 * num_pixels
    energy_indices = np.arange(energy_start, file_size, cluster_bytes)
    energy_bytes = np.stack([raw[energy_indices + j] for j in range(4 * num_pixels)], axis=0)
    energies = np.frombuffer(energy_bytes.T.tobytes(), dtype=np.float32)
    energies = energies.reshape((n_clusters, num_pixels)).T

    # --- Step 5. Return dictionary ---
    clusters = {
        'frame_ids': frame_ids,
        'pixel_ids': pixel_ids,
        'adus': adus,
        'energies': energies
    }

    t1 = time.time()
    print(f"Loaded {file_size/1e6:.2f} MB in {t1 - t0:.2f} s (≈ {file_size/1e6/(t1-t0+1e-9):.1f} MB/s)")

    return clusters

def write_clusters(filename, clusters):
    with open(filename, "wb") as file:
        for cluster in clusters:
            header = struct.pack("ii", cluster.frame_id, len(cluster.pixel_ids))
            file.write(header)
            file.write(struct.pack(f"{len(cluster.pixel_ids)}H", *cluster.pixel_ids))
            file.write(struct.pack(f"{len(cluster.adus)}H", *cluster.adus))
            file.write(struct.pack(f"{len(cluster.energies)}f", *cluster.energies))

if __name__ == "__main__":

    # # read and print clusters
    # shortened_clusters = read_clusters(r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_clusters\clusters_crystal_0_pixel_2.bin", n_max=3)
    # write_clusters(r"test_clusters_short.bin", shortened_clusters)
    # for i in range(len(shortened_clusters)):
    #     print(shortened_clusters[i])
    # print("Clusters written and read back successfully match.")

    # test fast reader
    # only applies to fixed-size clusters
    clusters_dict = read_clusters_fast(r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_clusters_backup\clusters_crystal_0_pixel_1.bin")
    for i in range(3):
        frame_id = clusters_dict['frame_ids'][i]
        pixel_ids = clusters_dict['pixel_ids'][:, i]
        adus = clusters_dict['adus'][:, i]
        energies = clusters_dict['energies'][:, i]
        cluster = Cluster(frame_id, pixel_ids, adus, energies)
        print(cluster)