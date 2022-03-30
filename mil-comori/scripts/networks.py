import torch

# MIL-CoMorI model specification

class VirtualCellClassifier(torch.nn.Module):
    def __init__(self,
                 n_input_features,
                 n_virtual_cells,
                 n_classes):
        super(VirtualCellClassifier,self).__init__()
        self.n_input_features = n_input_features
        self.n_virtual_cells = n_virtual_cells
        self.n_classes = n_classes
        self.f = torch.nn.Sequential(
            torch.nn.Conv2d(
                in_channels=1,
                out_channels=self.n_virtual_cells,
                kernel_size=[1,self.n_input_features],
                padding=0),
            torch.nn.Softmax(dim=1)
        )
        self.h = torch.nn.Sequential(
            torch.nn.Linear(self.n_virtual_cells,
                            self.n_classes),
            torch.nn.Softmax(dim=-1)
        )

    def g(self,x):
        x = x / x.shape[2]
        return torch.sum(x,dim=2).squeeze(dim=2)

    def forward(self,x):
        x = self.f(x)
        x = self.g(x)
        x = self.h(x)
        return x

    def get_virtual_cells(self,x):
        x = self.f(x)
        x = self.g(x)
        return x

class VirtualCellClassifierStack(torch.nn.Module):
    def __init__(self,
                 n_input_features,
                 n_virtual_cells,
                 other_datasets_sizes,
                 n_classes,
                 dropout_rate=0.):
        """
        Multi-objective stack of VCQ with other datasets.

        args:
            * n_input_features - list of number of input
            features
            * n_virtual_cells - list of number of virtual
            cells (len(n_virtual_cells) == len(n_input_features))
            * other_datasets_sizes - size of other datasets to be
            combined in the last layer
            * n_classes - list of number of classes
        """
        super(VirtualCellClassifierStack,self).__init__()
        self.n_input_features = n_input_features
        self.n_virtual_cells = n_virtual_cells
        self.other_datasets_sizes = other_datasets_sizes
        self.n_classes = n_classes
        self.dropout_rate = dropout_rate
        self.total_n_final = sum(n_virtual_cells) + sum(other_datasets_sizes)
        self.initialize_vcq()
        self.initialize_final_layers()

    def initialize_vcq(self):
        self.vcq_list = torch.nn.ModuleList()
        for nif,nvc in zip(self.n_input_features,self.n_virtual_cells):
            vcq = VirtualCellClassifier(nif,nvc,2)
            self.vcq_list.append(vcq)

    def initialize_final_layers(self):
        self.final_layers = torch.nn.ModuleList()
        if self.dropout_rate != 0:
            self.dropout = torch.nn.Dropout(self.dropout_rate)
        for n_class in self.n_classes:
            final_layer = torch.nn.Sequential(
                torch.nn.Linear(
                    self.total_n_final,
                    n_class),
                torch.nn.Softmax(dim=-1))
            self.final_layers.append(final_layer)

    def __getitem__(self,i):
        return self.vcq_list[i]

    def forward(self,x,i=0,dropout=False):
        """
        Inferrence for final_layer i
        args:
            * x - list with [[cells],[other_datasets]]. If no other
            datasets this should be [[cells],None].
            * i - index of the final_layer
            * dropout - whether to use dropout
        """
        cell_dataset,other_datasets = x
        features = []
        for vcq,cd in zip(self.vcq_list,cell_dataset):
            quantified_cells = vcq.get_virtual_cells(cd)
            if dropout == True:
                quantified_cells = self.dropout(quantified_cells)
            features.append(quantified_cells)
        if other_datasets is not None:
            features.extend(other_datasets)
        feature_tensor = torch.cat(features,dim=-1)
        return self.final_layers[i](feature_tensor)

    def forward_w_mask(self,x,i=0,masks=None):
        """
        Inferrence for final_layer i
        args:
            * x - list with [[cells],[other_datasets]]. If no other
            datasets this should be [[cells],None].
            * i - index of the final_layer
            * dropout - whether to use dropout
            * mask - mask for virtual cells
        """
        cell_dataset,other_datasets = x
        features = []
        for vcq,cd,mask in zip(self.vcq_list,cell_dataset,masks):
            quantified_cells = vcq.get_virtual_cells(cd) * mask
            quantified_cells /= quantified_cells.sum()
            features.append(quantified_cells)
        if other_datasets is not None:
            features.extend(other_datasets)
        feature_tensor = torch.cat(features,dim=-1)
        return self.final_layers[i](feature_tensor)

    def get_cells(self,x,i=0):
        """
        Virtual cell proportions for classifier i
        args:
            * x - Cell features
            * i - index of the final_layer
        """
        cell_dataset = x
        features = []
        for vcq,cd in zip(self.vcq_list,cell_dataset):
            features.append(vcq.get_virtual_cells(cd))
        return features

    def get_cells_individual(self,x):
        """
        Virtual cell classifications
        args:
            * x - Cell features
        """
        cell_dataset = x
        features = []
        for vcq,cd in zip(self.vcq_list,cell_dataset):
            features.append(vcq.f(cd))
        return features
