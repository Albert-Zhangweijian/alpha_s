import struct


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

def write_clusters(filename, clusters):
    with open(filename, "wb") as file:
        for cluster in clusters:
            header = struct.pack("ii", cluster.frame_id, len(cluster.pixel_ids))
            file.write(header)
            file.write(struct.pack(f"{len(cluster.pixel_ids)}H", *cluster.pixel_ids))
            file.write(struct.pack(f"{len(cluster.adus)}H", *cluster.adus))
            file.write(struct.pack(f"{len(cluster.energies)}f", *cluster.energies))

if __name__ == "__main__":
    clusters = read_clusters(r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\clusters_crystal_0.bin", n_max=10)
    shortened_clusters = clusters[:5]
    write_clusters(r"test_clusters_short.bin", shortened_clusters)
    reread_clusters = read_clusters(r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\clusters_crystal_0_Co57_reset_energies.bin", n_max=10)
    for i in range(len(shortened_clusters)):
        print(shortened_clusters[i])
        print(reread_clusters[i])
    print("Clusters written and read back successfully match.")
        