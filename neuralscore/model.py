"""EpiGuide NeuroScore model — SparseForcedEdgeGNN.

Self-contained architecture matching the released weights (model/model_state.pt).
A sparse graph neural network over a fixed CpG similarity graph with frozen,
task-supervised edge weights and masked pooling over observed CpG nodes.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing


class FixedEdgeWeightConv(MessagePassing):
    """Message passing where each edge's (frozen) scalar weight scales the message."""
    def __init__(self, hidden_dim, dropout=0.1):
        super().__init__(aggr="add")
        self.message_mlp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, hidden_dim),
        )
        self.update_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim), nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.norm = nn.LayerNorm(hidden_dim)

    def forward(self, x, edge_index, edge_attr):
        out = self.propagate(edge_index=edge_index, x=x, edge_attr=edge_attr)
        x_new = self.update_mlp(torch.cat([x, out], dim=1))
        return self.norm(x + x_new)

    def message(self, x_j, edge_attr):
        weight = edge_attr.view(-1, 1)
        return weight * self.message_mlp(x_j)


class SparseForcedEdgeGNN(nn.Module):
    """Node features (3): methylation call, observed mask, coverage. Output: logit for P(neural-high)."""
    def __init__(self, node_in_dim=3, hidden_dim=32, num_layers=2, dropout=0.1):
        super().__init__()
        self.node_encoder = nn.Linear(node_in_dim, hidden_dim)
        self.convs = nn.ModuleList(
            [FixedEdgeWeightConv(hidden_dim=hidden_dim, dropout=dropout) for _ in range(num_layers)]
        )
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, 1),
        )

    def forward(self, data):
        x = self.node_encoder(data.x)
        for conv in self.convs:
            x = conv(x, data.edge_index, data.edge_attr)
            x = F.relu(x)

        observed_mask = data.x[:, 1:2]
        batch = data.batch
        n = batch.max().item() + 1

        # masked mean over observed nodes
        x_obs = x * observed_mask
        pooled_obs = torch.zeros(n, x.size(1), device=x.device).index_add(0, batch, x_obs)
        counts = torch.zeros(n, 1, device=x.device).index_add(0, batch, observed_mask)
        pooled_obs = pooled_obs / counts.clamp(min=1.0)

        # mean over all nodes
        pooled_all = torch.zeros(n, x.size(1), device=x.device).index_add(0, batch, x)
        node_counts = torch.bincount(batch).float().view(-1, 1).to(x.device)
        pooled_all = pooled_all / node_counts.clamp(min=1.0)

        pooled = torch.cat([pooled_obs, pooled_all], dim=1)
        return self.classifier(pooled).squeeze(-1)
