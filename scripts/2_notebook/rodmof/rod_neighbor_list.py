import pormake

class RodNeighborList(pormake.neighbor_list.NeighborList):
    def __init__(self, n_list):
        self._neighbor_list = n_list
