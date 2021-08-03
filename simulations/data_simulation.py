import numpy as np
from sklearn.datasets import make_spd_matrix
import torch

class DataGenerator(torch.utils.data.Dataset):
    def __init__(self,
                 n_features,n_virtual_cells,
                 n_classes,n_cells,batch_size,
                 standard_deviation_centers,
                 shared_variance=True,
                 multivariate_normal=False,
                 n_cpus=1):
        self.n_features = n_features
        self.n_virtual_cells = n_virtual_cells
        self.n_classes = n_classes
        self.n_cells = n_cells
        self.batch_size = batch_size
        self.standard_deviation_centers = standard_deviation_centers
        self.shared_variance = shared_variance
        self.multivariate_normal = multivariate_normal
        self.n_cpus = n_cpus

        self.generate_initial_parameters()

    def generate_initial_parameters(self):
        self.cell_class_centers = np.random.normal(
            0,10,
            size=[self.n_virtual_cells,self.n_features])

        if self.shared_variance == False:
            self.variance = np.random.exponential(size=[self.n_virtual_cells,
                                                        self.n_features])
        else:
            self.variance = np.ones(
                [self.n_virtual_cells,self.n_features]) * np.random.exponential(size=1)

        if self.multivariate_normal == True:
            self.random_cov_dist = np.stack(
                [make_spd_matrix(self.n_features)
                 for _ in range(self.n_virtual_cells)],
                axis=0)

            self.mvn_dist = [torch.distributions.MultivariateNormal(
                torch.Tensor(self.cell_class_centers[i,:]),
                torch.Tensor(self.random_cov_dist[i,:,:]))
                for i in range(self.n_virtual_cells)]
        else:
            self.n_dist = [torch.distributions.Normal(
                torch.Tensor(self.cell_class_centers[i,:]),
                torch.Tensor(self.variance[i,:]))
                for i in range(self.n_virtual_cells)]

        self.class_proportions = [
            np.random.dirichlet(
                np.ones(self.n_virtual_cells)/10) for _ in range(self.n_classes)
        ]

    def generate_item(self,curr_class):
        member_distribution = np.random.multinomial(
            self.n_cells,
            pvals=self.class_proportions[curr_class]
        )
        population = []
        for i,m in enumerate(member_distribution):
            if self.multivariate_normal == False:
                P = self.n_dist[i].sample((m,))
            else:
                P = self.mvn_dist[i].sample((m,))
            population.append(P)
        population = np.concatenate(population,axis=0)
        return population

    def generate_batch(self,batch_size=None):
        if not batch_size:
            batch_size = self.batch_size
        batch = []
        truth = []
        for _ in range(batch_size):
            curr_class = np.random.randint(len(self.class_proportions))
            population = self.generate_item(curr_class)
            batch.append(population)
            truth.append(curr_class)
        batch = np.stack(batch,axis=0)
        truth = np.array(truth)
        return batch,truth

    def __call__(self):
        return self.generate_batch()

    def __getitem__(self,i):
        """
        This method is implemented for the sole purpose of using
        PyTorch dataloaders.
        """
        i = i % self.n_classes
        item = self.generate_item(i)
        return item,i

    def __len__(self):
        """
        This method is implemented for the sole purpose of using
        PyTorch dataloaders.
        """
        return int(1e6 * self.n_classes)
