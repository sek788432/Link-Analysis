import numpy as np
class link_analysis:
    def __init__(self, node_num=0, adj_mtx=None, node_dict={}, page_rk=None, d_factor = 0.15, hub=None, authority=None, epsilon = 1e-6, sim_rk = None):
        self.node_num = node_num
        self.adj_mtx = adj_mtx
        self.node_dict = node_dict
        self.page_rk = page_rk
        self.d_factor = d_factor
        self.hub = hub
        self.authority = authority
        self.epsilon = epsilon
        self.sim_rk = sim_rk

    def get_map(self, file_name):
        file = open(file_name)
        data = file.readlines()
        file.close()
        nodes = set()
        edge = []
        for line in data:
            f_node, s_node = line.split(",")
            s_node = s_node.rstrip()
            edge.append([f_node,s_node])
            nodes.add(f_node)
            nodes.add(s_node)
        self.node_num = len(nodes)
        self.page_rk = [1/self.node_num for i in range(self.node_num)]
        self.hub = [1 for i in range(self.node_num)]
        self.authority = [1 for i in range(self.node_num)]
        nodes = sorted(nodes)
        for i in range(len(nodes)):
            self.node_dict[nodes[i]] = i
        self.adj_mtx = np.zeros((self.node_num, self.node_num))
        for i in edge:
            self.adj_mtx[self.node_dict[i[0]], self.node_dict[i[1]]] = 1

    def cal_pagerank(self):
        while True:
            tmp_page_rk = []
            for i in range(self.node_num):
                tmp_page_rk.append(self.d_factor/self.node_num)     #first half
                sum = 0
                for j in range(self.node_num):
                    if self.adj_mtx[j][i] == 1:
                        out_link = 0
                        for k in range(self.node_num):
                            if self.adj_mtx[j][k] == 1:
                                out_link += 1
                        sum += self.page_rk[j]/out_link
                tmp_page_rk[i] += sum*(1-self.d_factor)             #second half
            if np.linalg.norm(np.asarray(tmp_page_rk) - np.asarray(self.page_rk) , 2) < self.epsilon:
                break
            self.page_rk = tmp_page_rk.copy()

    def HITS(self):
        while True:
            authority = np.dot(self.adj_mtx.T, np.asarray(self.hub))
            hub = np.dot(self.adj_mtx, authority)
            authority /= np.linalg.norm(authority,1)    # l1 normalize
            hub /= np.linalg.norm(hub,1)                # l1 normalize
            a_dist = np.linalg.norm(authority - self.authority, 2)
            h_dist = np.linalg.norm(hub - self.hub, 2)
            self.hub = hub
            self.authority = authority
            if (a_dist + h_dist) < self.epsilon:
                break

    def simrank(self, c, iteration):

        graph = self.adj_mtx
        row, col = graph.shape
        identity = np.identity(row)   # Init SimRank
        parents = []                  # Parents of each node
        num_p = []                    # Num of parents of each node
        for i in range(row):
            tmp = np.nonzero(graph[:, i])[0]
            parents.append(tmp)
            num_p.append(tmp.shape[0])

        for k in range(iteration):
            new_simrk = np.zeros((row, col))   # New SimRank for new iteration
            for a in range(row):
                for b in range(col):
                    if a == b:    # s(a, a) = 1
                        new_simrk[a, b] = 1
                        continue
                    if num_p[a] == 0 or num_p[b] == 0:  # s(a, b) = 0 if a or b has no parent
                        new_simrk[a, b] = 0
                        continue
                    tmp = 0
                    for i in range(num_p[a]):
                        for j in range(num_p[b]):
                            tmp += identity[parents[a][i], parents[b][j]]
                    new_simrk[a, b] = c / (num_p[a] * num_p[b]) * tmp
            self.sim_rk = new_simrk

    def make_txt(self, file_name):
        folder_path = "./output_file/"
        np.savetxt(folder_path + file_name + "_HITS_hub.txt", self.hub, newline = " ", fmt='%f')
        np.savetxt(folder_path + file_name + "_HITS_authority.txt", self.authority, newline = " ", fmt='%f')
        np.savetxt(folder_path + file_name + "_PageRank.txt", self.page_rk, newline = " ", fmt='%f')
        np.savetxt(folder_path + file_name + "_SimRank.txt", self.sim_rk, fmt='%f')


if __name__ == '__main__':
    file_name = "./dataset/graph_1.txt"
    graph = link_analysis()
    graph.get_map(file_name)
    graph.cal_pagerank()
    graph.HITS()
    graph.simrank(0.6, 10)
    graph.make_txt(file_name[13:-4])
    print("Hub:\n", graph.hub)
    print("Authority:\n", graph.authority)
    print("PageRank:\n", graph.page_rk)
    print("SimRank:\n", graph.sim_rk)
