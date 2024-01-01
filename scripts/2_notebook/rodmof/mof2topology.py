import os
import copy
import math
from pathlib import Path
from collections import defaultdict

import numpy as np
import networkx as nx

import ase
import ase.visualize
import ase.neighborlist

import pormake

from .rod_utils import *
from .rod_neighbor_list import RodNeighborList



class MOF2TOPO:
    def __init__(self, atoms):
        self.atoms = atoms
        self.cleanup()
        self.find_building_blocks()
        self.bb_found = False
        self.coord_check = []

        self._node_origins = dict()
        self._neighbor_nodes = dict()
        for node in self.nodes:
            node_origin = self.get_bb_center(node)
            self._node_origins[self._node_index(node)] = self.atoms.cell.scaled_positions(node_origin)
            self._neighbor_nodes[self._node_index(node)] = self.get_neighbor_nodes(node)



    def cleanup(self, remove_interpenetration=False):
        # Get bond except metals.
        I, J, _ = covalent_neigbor_list(self.atoms)
        # Build MOF graph.
        graph = nx.Graph(zip(I, J))

        ccs = sorted(nx.connected_components(graph), reverse=True, key=len)

        if len(ccs) < 2:
            indices = list(range(len(self.atoms)))
        elif remove_interpenetration:
            indices = list(ccs[0])
        elif len(ccs[0]) == len(ccs[1]):
            indices = list(ccs[0] | ccs[1])
        else:
            indices = list(ccs[0])
        self.atoms = self.atoms[indices]


    def find_building_blocks(self):
        # Get full bond information.
        I, J, _ = covalent_neigbor_list(self.atoms)
        bond_list = [[] for _ in range(len(self.atoms))]
        # Build neighbor list as a list form.
        for i, j in zip(I, J):
            bond_list[i].append(j)

        # Get indices of metal atoms.
        metal_indices = [i for i, a in enumerate(self.atoms) if a.symbol in METAL_LIKE]

        # Mark liking atom indices.
        liking_atom_indices = []
        for i in range(len(self.atoms)):
            if set(bond_list[i]) & set(metal_indices):
                liking_atom_indices.append(i)
        
        # Build MOF graph.
        graph = nx.Graph(zip(I, J))
        # Remove metal containing edges.
        metal_containing_edges = list(graph.edges(metal_indices))

        test_graph = graph.copy()
        test_graph.remove_edges_from(metal_containing_edges)
        result = []
        for cc in list(nx.connected_components(test_graph)):
            # Neglect single node components.
            if len(cc) == 1:
                continue

            # Construct graph of connected component.
            cc_graph = nx.subgraph(graph, cc).copy()

            # Get all bridges.
            bridges = list(nx.bridges(cc_graph))

            # Filter bridges.
            # This filter not filter out the self liking bridges.
            filtered_bridges = []
            for b in bridges:
                test_graph = cc_graph.copy()
                test_graph.remove_edge(*b)
                c1, c2 = list(nx.connected_components(test_graph))

                # Neglect no metal components.
                if not set(liking_atom_indices) & c1:
                    continue
                elif not set(liking_atom_indices) & c2:
                    continue

                # Neglect if it is not connected to metal
                elif len(c1) == 1 and (c1 not in liking_atom_indices):
                    continue
                elif len(c2) == 1 and (c2 not in liking_atom_indices):
                    continue

                filtered_bridges.append(b)

            # Get first level bridges only.
            test_graph = cc_graph.copy()
            test_graph.remove_edges_from(filtered_bridges)
            test_ccs = list(nx.connected_components(test_graph))

            liking_ccs = []
            for test_cc in test_ccs:
                if set(liking_atom_indices) & test_cc:
                    liking_ccs.append(test_cc)

            first_level_bridges = set()
            for liking_cc in liking_ccs:
                for b in filtered_bridges:
                    if set(b) & liking_cc:
                        first_level_bridges.add(b)
            first_level_bridges = list(first_level_bridges)

            result += first_level_bridges
        first_level_bridges = result

        # Remove self liking ligands.
        test_graph = graph.copy()
        #self.connecting_site_list = np.unique(first_level_bridges)
        test_graph.remove_edges_from(first_level_bridges)
        building_blocks = list(nx.connected_components(test_graph))
        #self._building_blocks = building_blocks

        # Merge self connecting linkers that form a path of bb to the same bb.
        index_to_bb = {}
        for i, bb in enumerate(building_blocks):
            for j in bb:
                index_to_bb[j] = i

        merging_dict = collections.defaultdict(list)
        for i, bb in enumerate(building_blocks):
            species = set(self.atoms[list(bb)].symbols)
            if species & set(METAL_LIKE):
                continue
            # Get connection point.
            connection_indices = []
            for j in bb:
                for k in graph.adj[j].keys():
                    if index_to_bb[k] == i:
                        continue
                    connection_indices.append(k)
            linked_bb_indices = [index_to_bb[_] for _ in connection_indices]
            if len(set(linked_bb_indices)) == 1:
                parent_bb_index = linked_bb_indices[0]
                merging_dict[parent_bb_index].append(i)

        children_indices = []
        new_bb = copy.deepcopy(building_blocks)
        for k, v in merging_dict.items():
            # Save index of child bb to remove later.
            tobemerged = set()
            for bb_index in v:
                tobemerged |= building_blocks[bb_index]
            new_bb[k] |= tobemerged
            # Check dimension changes.
            original_atoms = self.atoms[list(building_blocks[k])]
            original_dim = evaluate_atoms_dimension(original_atoms)

            new_atoms = self.atoms[list(new_bb[k])]
            new_dim = evaluate_atoms_dimension(new_atoms)

            if original_dim == new_dim:
                # Accept merging.
                children_indices += v
            else:
                # Reject merging.
                new_bb[k] = building_blocks[k]

        new_bb = [bb for i, bb in enumerate(new_bb) if i not in children_indices]

        # Update building_block.
        building_blocks = new_bb

        # Remove improper first level bridges.
        index_to_bb = {}
        for i, bb in enumerate(building_blocks):
            for j in bb:
                index_to_bb[j] = i

        def check(bridge):
            i, j = bridge
            return index_to_bb[i] != index_to_bb[j]

        first_level_bridges = [b for b in first_level_bridges if check(b)]

        self._connecting_bonds = first_level_bridges
        self._connecting_site_list = np.unique(first_level_bridges)
        self._building_blocks = new_bb
        self.bb_found = True

        self._nodes = []
        self._edges = []
        for bb in self._building_blocks:
            connection_count = len(set(self._connecting_site_list) & bb)
            if connection_count == 1:
                raise InputError("[ERROR] node connection == 1")
            elif connection_count == 2:
                self._edges.append(bb)
            else:
                self._nodes.append(bb)


    @property
    def building_blocks(self):
        if not self.bb_found:
            self.find_building_blocks()
        return self._building_blocks


    @property
    def nodes(self):
        if not self.bb_found:
            self.find_building_blocks()
        return self._nodes


    @property
    def edges(self):
        if not self.bb_found:
            self.find_building_blocks()
        return self._edges


    @property
    def connecting_site_list(self):
        if not self.bb_found:
            self.find_building_blocks()
        return self._connecting_site_list


    @property
    def connecting_bonds(self):
        if not self.bb_found:
            self.find_building_blocks()
        return self._connecting_bonds


    def get_bb_center(self, bb):
        # from rod_utils
        bb_atoms = self.atoms[list(bb)]
        return get_centroid(bb_atoms)


    def get_neighbor_bbs(self, bb):

        neighbors = []
        # connection_indices : inside bb
        # target_indices : outside bb
        connection_indices = set(self._connecting_site_list) & bb
        target_indices = []
        for i in connection_indices:
            for bond in self.connecting_bonds:
                if i in bond:
                    if bond[0] == i:
                        target_indices.append(bond[1])
                    elif bond[1] == i:
                        target_indices.append(bond[0])
                    break

        target_indices = set(target_indices)
        
        # len(neighbors) == len(target_indices)
        def unit_vector(vector):
            """ Returns the unit vector of the vector.  """
            return vector / np.linalg.norm(vector)

        def angle_between(v1, v2):
            """ Returns the angle in radians between vectors 'v1' and 'v2' """
            v1_u = unit_vector(v1)
            v2_u = unit_vector(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


        # CASE : EDGE
        if bb in self.edges:
            # Rod-part (need to fix?)
            if len(bb) == 1:
                for b in self.nodes:
                    if b != bb:
                        if target_indices & b:
                            n_center = self.get_bb_center(b)
                            center = self.get_bb_center(bb)
                            n_vector = n_center - center
                            n_vector = self.atoms.cell.scaled_positions(n_vector)
                            for i in range(3):
                                if n_vector[i] > 0.5:
                                    n_vector[i] -= 1.0
                                if n_vector[i] < -0.5:
                                    n_vector[i] += 1.0
                            target_bond = [list(bb)[0], list(bb)[0]]
                            neighbors.append([b, n_vector, target_bond])

                neighbor_sum = np.sum([j for _, j, _ in neighbors], axis=0)
                for i in range(3):
                    if neighbor_sum[i] > 0.75:
                        cand = sorted(neighbors, key=lambda x: np.linalg.norm(x[1]))[-1][0]
                        for c, vec in neighbors:
                            if cand == c:
                                vec[i] -= 1.0
                    if neighbor_sum[i] < -0.75:
                        cand = sorted(neighbors, key=lambda x: np.linalg.norm(x[1]))[-1][0]
                        for c, vec in neighbors:
                            if cand == c:
                                vec[i] += 1.0
                return neighbors


            for index in target_indices:
                # find target bb
                target = [b for b in self.nodes if index in b][0]
                # edge_center
                edge_center = self.get_bb_center(bb)
                # direction_vector : target - edge center
                direction_vec = self.atoms.get_positions()[index] - edge_center
                direction_vec = self.atoms.cell.scaled_positions(direction_vec)
                for i in range(3):
                    if direction_vec[i] > 0.5:
                        direction_vec[i] -= 1.0
                    if direction_vec[i] < -0.5:
                        direction_vec[i] += 1.0

                # target_center
                target_center = self.get_bb_center(target)
                target_vec = target_center - edge_center
                target_vec = self.atoms.cell.scaled_positions(target_vec)
                for i in range(3):
                    if target_vec[i] > 0.5:
                        target_vec[i] -= 1.0
                    if target_vec[i] < -0.5:
                        target_vec[i] += 1.0

                min_angle = 180
                min_vec = None
                for i in [0,1,-1]:
                    for j in [0,1,-1]:
                        for k in [0,1,-1]:
                            vec = target_vec + np.array([i,j,k])
                            if any(vec > 0.7):
                                continue
                            if any(vec < -0.7):
                                continue
                            angle = angle_between(vec, direction_vec)
                            if angle < min_angle:
                                min_angle = angle
                                min_vec = vec

                target_bond = [bond for bond in self.connecting_bonds if index in bond][0]

                neighbors.append([target, min_vec, target_bond])

        elif bb in self.nodes:
            for index in target_indices:
                # find target bb
                target = [b for b in self.edges if index in b][0]
                target_bond = [bond for bond in self.connecting_bonds if index in bond][0]
                target_neighbors = self.get_neighbor_bbs(target)
                for n, vec, bond in target_neighbors:
                    if all([tuple(bond) == target_bond, bb == n]):
                        neighbors.append([target, -vec, bond])

        else:
            print("BB Error")
        return neighbors


    def build_topology_atoms(self):
        seqs = self.get_coord_seq()
        seq_map = dict()
        for i, seq in enumerate(seqs):
            seq_map[i] = list(seq)[:7]

        spacegroup = 'P1'
        node_positions = []
        coordination_numbers = []
        node_tags = []
        for i, node in enumerate(self.nodes):
            seq = self.get_coord_seq_from_node(node, length=7, test_length=11)
            for j in seq_map.keys():
                if seq_map[j] == seq:
                    node_tags.append(j)
            pos = self.get_bb_center(node)
            pos = self.atoms.cell.scaled_positions(pos)
            for i in range(3):
                if pos[i] > 0.5:
                    pos[i] -= 1.0
                if pos[i] < -0.5:
                    pos[i] += 1.0

            node_positions.append(pos)
            connection_count = len(set(self._connecting_site_list) & node)
            coordination_numbers.append(connection_count)

        node_positions = np.array(node_positions)

        edge_center_positions = []
        for edge in self.edges:
            pos = self.get_bb_center(edge)
            pos = self.atoms.cell.scaled_positions(pos)
            for i in range(3):
                if pos[i] > 0.5:
                    pos[i] -= 1.0
                if pos[i] < -0.5:
                    pos[i] += 1.0

            edge_center_positions.append(pos)
        edge_center_positions = np.array(edge_center_positions)

        n_nodes = node_positions.shape[0]
        n_edges = edge_center_positions.shape[0]
        species = np.concatenate([
            np.full(shape=n_nodes, fill_value="C"),
            np.full(shape=n_edges, fill_value="O"),
        ])
        coords = np.concatenate([node_positions, edge_center_positions], axis=0)

        node_types = node_tags #[i for i,_ in enumerate(node_positions)]
        edge_types = [-1 for _ in edge_center_positions]
        site_properties = {
            "type": node_types + edge_types,
            "cn": coordination_numbers + [2 for _ in edge_center_positions],
        }

        structure = mg.Structure.from_spacegroup(
            sg='P1',
            lattice=mg.Lattice.from_parameters(*self.atoms.get_cell_lengths_and_angles()),
            species=species,
            coords=coords,
            site_properties=site_properties,
        )

        atoms = ase.Atoms(
            symbols=[s.name for s in structure.species],
            positions=structure.cart_coords,
            cell=structure.lattice.matrix,
            tags=structure.site_properties["type"],
            pbc=True,
        )
        return atoms


    def _periodic_count(self, bb, pos):
        bb_origin = self.get_bb_center(bb)
        bb_origin = self.atoms.cell.scaled_positions(bb_origin)
        return [int(i) for i in np.round(pos - bb_origin)]





    def get_neighbor_nodes(self, bb):
        n_list = []
        for i, v1, b1 in self.get_neighbor_bbs(bb):
            for j, v2, b2 in self.get_neighbor_bbs(i):
                if not all([j==bb, b1==b2]):
                    vec = v1 + v2
                    n_list.append([j, vec])
                    break
        return n_list


    def _node_index(self, bb):
        for i, b in enumerate(self.nodes):
            if b == bb:
                return i


    def _bb_index(self, bb):
        for i, b in enumerate(self.building_blocks):
            if b == bb:
                return i


    def get_coord_seq_from_node(self, node, length=10, test_length=7):
        seq = []
        used = dict()
        for i in range(len(self.nodes)):
            used[i] = []
        end = [(node, [0,0,0])]
        used[self._node_index(node)].append([0,0,0])

        for i in range(length):
            new_end = []
            for n, index in end:
                node_origin = self._node_origins[self._node_index(n)]
                pos = node_origin + np.array(index)
                for b, vector in self._neighbor_nodes[self._node_index(n)]:
                    b_pos = pos + vector
                    b_index = self._periodic_count(b, b_pos)
                    if not b_index in used[self._node_index(b)]:
                        new_end.append([b, b_index])
                        used[self._node_index(b)].append(b_index)

            end = copy.deepcopy(new_end)
            seq.append(len(end))

            if len(seq) == test_length:
                if seq in self.coord_check:
                    return []
                else:
                    self.coord_check.append(copy.deepcopy(seq))
        return seq


    def get_coord_seq(self):
        seqs = []
        for node in self.nodes:
            seq = self.get_coord_seq_from_node(node)
            if seq:
                #print('check : ', seq)
                seqs.append(seq)

        return np.array(sorted(seqs))


    def get_pormake_neighbor_list(self):
        n_list = [[] for _ in range(len(self.building_blocks))]
        for bb in self.building_blocks:
            neighbors = self.get_neighbor_bbs(bb)
            bb_index = self._bb_index(bb)
            for b, vec, _ in neighbors:
                b_index = self._bb_index(b)
                n_list[bb_index].append(pormake.neighbor_list.Neighbor(b_index, vec))

        return RodNeighborList(n_list)
