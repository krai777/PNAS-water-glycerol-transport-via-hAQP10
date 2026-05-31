import MDAnalysis as mda
import numpy as np
from collections import defaultdict


class PermeationTracker:
    def __init__(self, psf_file, dcd_file, cylinder_radius=8.0, z_min=-18.0, z_max=18.0):
        self.u = mda.Universe(psf_file, dcd_file)
        self.cylinder_radius = cylinder_radius
        self.z_min = z_min
        self.z_max = z_max

        # Atom index ranges for each monomer chain
        self.chain_indices = {
            'A': (282, 4040),
            'B': (4802, 8574),
            'C': (9278, 13080),
            'D': (13739, 17600)
        }

        # chain atom selections once
        self.chain_atoms = {
            chain: self.u.select_atoms(f"index {start_idx}:{end_idx}")
            for chain, (start_idx, end_idx) in self.chain_indices.items()
        }

        # Select permeating species
        # For glycerol example:
        # self.waters = self.u.select_atoms("resname MGLY and name C1")
        self.waters = self.u.select_atoms("resname TIP3 and name OH2")

        self.chain_data = {
            chain: {
                'water_tracking': defaultdict(lambda: {
                    'inside': False,
                    'entry_side': None,
                    'entry_frame': None,
                    'crossed_center': False,
                    'fully_exited': False,
                    'last_position': None,   # now stores relative z
                    'permeation_count': 0,
                    'entry_exit_history': []
                }),
                'permeation_counts': {'up': 0, 'down': 0},
                'frame_counts': [],
                'permeation_details': []
            } for chain in self.chain_indices.keys()
        }

    def is_in_cylinder(self, point, com):
        """
        Check whether a point lies inside the cylindrical pore region
        defined relative to the chain COM.
        """
        dx = point[0] - com[0]
        dy = point[1] - com[1]
        dz = point[2] - com[2]

        dist_xy = np.sqrt(dx * dx + dy * dy)
        return (dist_xy <= self.cylinder_radius) and (self.z_min <= dz <= self.z_max)

    def get_side(self, rel_z):
        """
        Determine whether the molecule is on the upper or lower side
        relative to the chain COM plane.
        """
        return 'upper' if rel_z > 0 else 'lower'

    def analyze_frame(self, ts):
        for chain in self.chain_indices.keys():
            protein_chain = self.chain_atoms[chain]
            com = protein_chain.center_of_mass()

            frame_up = 0
            frame_down = 0

            for water in self.waters:
                water_id = water.index
                water_pos = water.position
                rel_z = water_pos[2] - com[2]

                tracking = self.chain_data[chain]['water_tracking'][water_id]

                if self.is_in_cylinder(water_pos, com):
                    current_side = self.get_side(rel_z)

                    # Molecule enters the pore region
                    if not tracking['inside']:
                        tracking['inside'] = True
                        tracking['entry_side'] = current_side
                        tracking['entry_frame'] = ts.frame
                        tracking['crossed_center'] = False
                        tracking['fully_exited'] = False
                        tracking['last_position'] = rel_z

                    # Check whether molecule crossed the central plane
                    if tracking['last_position'] is not None:
                        if (tracking['last_position'] < 0 and rel_z > 0) or \
                           (tracking['last_position'] > 0 and rel_z < 0):
                            tracking['crossed_center'] = True

                    tracking['last_position'] = rel_z

                else:
                    # Molecule is outside pore region
                    if tracking['inside'] and tracking['crossed_center'] and not tracking['fully_exited']:
                        exit_side = self.get_side(rel_z)

                        # Count only if it exits on the opposite side beyond z-limits
                        if tracking['entry_side'] != exit_side and (rel_z < self.z_min or rel_z > self.z_max):
                            direction = 'up' if tracking['entry_side'] == 'lower' else 'down'

                            if direction == 'up':
                                frame_up += 1
                            else:
                                frame_down += 1

                            self.chain_data[chain]['permeation_details'].append({
                                'water_id': water_id,
                                'entry_frame': tracking['entry_frame'],
                                'exit_frame': ts.frame
                            })

                            tracking['permeation_count'] += 1
                            tracking['fully_exited'] = True

                    # Fully reset state only when molecule is clearly outside z-bounds
                    if rel_z < self.z_min or rel_z > self.z_max:
                        tracking['inside'] = False
                        tracking['entry_side'] = None
                        tracking['entry_frame'] = None
                        tracking['crossed_center'] = False
                        tracking['last_position'] = None
                        tracking['fully_exited'] = False

            self.chain_data[chain]['permeation_counts']['up'] += frame_up
            self.chain_data[chain]['permeation_counts']['down'] += frame_down

            self.chain_data[chain]['frame_counts'].append({
                'frame': ts.frame,
                'up_to_down': frame_up,
                'down_to_up': frame_down,
                'total_permeation': frame_up + frame_down,
                'cumsum_permeation': (
                    self.chain_data[chain]['permeation_counts']['up'] +
                    self.chain_data[chain]['permeation_counts']['down']
                )
            })

    def run_analysis(self, start=None, stop=None, step=None):
        for ts in self.u.trajectory[start:stop:step]:
            if ts.frame % 100 == 0:
                print(f"Processing frame {ts.frame}")
            self.analyze_frame(ts)

        self.save_results()

    def save_results(self):
        for chain in self.chain_indices.keys():
            with open(f'NEWwatcounts_{chain}.txt', 'w') as f:
                f.write("FRAME\tup_to_down\tdown_to_up\ttotal_permeation\tcumsum_permeation\n")
                for count in self.chain_data[chain]['frame_counts']:
                    f.write(
                        f"{count['frame']}\t"
                        f"{count['up_to_down']}\t"
                        f"{count['down_to_up']}\t"
                        f"{count['total_permeation']}\t"
                        f"{count['cumsum_permeation']}\n"
                    )

            with open(f'NEWwatdetails_{chain}.txt', 'w') as f:
                f.write("MoleculeID\tentry_frame\texit_frame\n")
                for event in self.chain_data[chain]['permeation_details']:
                    f.write(
                        f"{event['water_id']}\t"
                        f"{event['entry_frame']}\t"
                        f"{event['exit_frame']}\n"
                    )

    def get_results(self):
        return self.chain_data


if __name__ == "__main__":
    psf_file = "../../hAQP10_wi.psf" #psf of your system
    dcd_file = "../md.dcd" #wrapped trajectory file

    tracker = PermeationTracker(
        psf_file=psf_file,
        dcd_file=dcd_file,
        cylinder_radius=8.0,
        z_min=-18.0,
        z_max=18.0
    )

    tracker.run_analysis()
    results = tracker.get_results()

    print("\nPermeation Summary:")
    for chain in results.keys():
        up = results[chain]['permeation_counts']['up']
        down = results[chain]['permeation_counts']['down']
        total = up + down
        print(f"Chain {chain}: {up} up, {down} down, {total} total permeations")
