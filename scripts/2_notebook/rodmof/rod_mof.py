from .mof2topology import *
from pormake import Topology
import pormake
from pormake import *

class RodMOF(MOF2TOPO):
    def __init__(self, atoms, node_BBs, edge_BBs, bonds):
        self.atoms = atoms
        self.cleanup()
        self.node_BBs, self.edge_BBs = node_BBs, edge_BBs
        self.bb_found = True
        self._building_blocks = node_BBs + edge_BBs
        self._nodes = node_BBs
        self._edges = edge_BBs

        self._connecting_bonds = bonds
        self._connecting_site_list = np.unique(bonds)
        self.coord_check = []

        self._node_origins = dict()
        self._neighbor_nodes = dict()
        for node in self.nodes:
            node_origin = self.get_bb_center(node)
            self._node_origins[self._node_index(node)] = self.atoms.cell.scaled_positions(node_origin)
            self._neighbor_nodes[self._node_index(node)] = self.get_neighbor_nodes(node)


class RodTopoView(Topology):
    def __init__(self, atoms, name, n_list, local_func=None):
        self.atoms = atoms
        self.name = name
        self.neighbor_list = n_list
        self.calculate_properties()
        self.local_structure_func = local_func


class RodTopology(Topology):
    def add_rod_info(self, indices):
        self.rod_edge_indices = indices

    def get_rod_permutation_info(self):
        """
        {n : [T,F,T,F]}     n -> slot, 
                            True -> Rod,
                            False-> Edge,
        """
        topo_permutations = dict()
        for i in self.node_indices:
            rod_check = [k in self.rod_edge_indices for k in self.get_neighbor_indices(i)]
            topo_permutations[i] = rod_check
        return topo_permutations


class RodBuildingBlock(BuildingBlock):
    def add_rod_info(self, indices):
        self.rod_indices = indices

    def get_permutation_for_slot(self, topo_info):
        """
        topo_info = list
        return [0, .., .., n]
        """
        permutation = [None] * len(topo_info)
        rod_count = sum(topo_info)
        
        bb_rod_indices = []
        bb_edge_indices = []
        for i, k in enumerate(self.connection_point_indices):
            if k in self.rod_indices:
                bb_rod_indices.append(i)
            else:
                bb_edge_indices.append(i)

        bb_rod_indices.reverse()
        bb_edge_indices.reverse()

        for i, t in enumerate(topo_info):
            if t:
                # rod
                index = bb_rod_indices.pop()
                permutation[i] = index
            else:
                # edge
                index = bb_edge_indices.pop()
                permutation[i] = index

        return permutation
