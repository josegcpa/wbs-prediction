import torch

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
            #torch.nn.Sigmoid(),
        )
        self.g = lambda x: torch.mean(x,dim=2).squeeze(dim=2)
        self.h = torch.nn.Sequential(
            torch.nn.Linear(self.n_virtual_cells,
                            self.n_classes),
            torch.nn.Softmax(dim=-1)
        )

    def forward(self,x):
        x = self.f(x)
        x = self.g(x)
        x = self.h(x)
        return x

    def get_virtual_cells(self,x):
        return self.g(self.f(x))
