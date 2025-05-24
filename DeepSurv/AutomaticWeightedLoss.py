#!/usr/bin/env python
# coding=utf-8
'''
Author       : Bowen Zheng
Date         : 2023-07-10 11:27:52
LastEditors  : ibowennn shihun44@163.com
LastEditTime : 2023-07-21 21:04:51
Description  : 

'''
# -*- coding: utf-8 -*-

import torch
import torch.nn as nn
import torch.nn.functional as F


class AutomaticWeightedLoss(nn.Module):
    """automatically weighted multi-task loss

    Params：
        num: int，the number of loss
        x: multi-task loss
    Examples：
        loss1=1
        loss2=2
        awl = AutomaticWeightedLoss(2)
        loss_sum = awl(loss1, loss2)
    """

    def __init__(self, num=2):
        super(AutomaticWeightedLoss, self).__init__()
        params = torch.ones(num, requires_grad=True)
        self.params = torch.nn.Parameter(params)

    def forward(self, *x):
        loss_sum = 0
        for i, loss in enumerate(x):
            loss_sum += 0.5 / (self.params[i] ** 2) * loss + torch.log(
                1 + self.params[i] ** 2
            )
            # loss_sum += loss / (2 * self.params[i] ** 2) + torch.log(self.params[i])
        return loss_sum


def compute_grad_l2_norm(layers, loss):
    G = torch.autograd.grad(loss, layers, retain_graph=True, create_graph=True)
    G_norm = torch.cat([torch.norm(g, 2).unsqueeze(0) for g in G]).sum()
    return G_norm


class SimpleGradNormalizer:
    def __init__(self, lr_init=0.025, alpha=0.16, task_num=2):
        self.alpha = alpha
        self.init_loss = None
        # self.loss_weight = nn.Parameter(torch.ones(3, device='cuda:0'))
        # self.optim = torch.optim.Adam([self.loss_weight], lr=lr_init)
        default_device = 'cuda' if torch.cuda.is_available() else 'cpu'
        # default_device = 'cpu'
        self.device = default_device
        self.loss_weight = {
            f"task_{i}": nn.Parameter(torch.tensor([1.0], device=default_device))
            for i in range(task_num)
        }
        self.optim = torch.optim.Adam(
            [self.loss_weight[k] for k, _ in self.loss_weight.items()], lr=lr_init
        )

    def set_init_loss(self, losses):
        self.init_loss = {
            n: torch.tensor([l.item()], device=l.device) for n, l in losses.items()
        }

    def normalize_loss_weight(self):
        num_losses = len(self.init_loss)
        # self.loss_weight.data[:] = torch.clamp(self.loss_weight, min=0.0)
        [
            self.loss_weight[n].data[:].clamp_(min=0.0)
            for n, _ in self.loss_weight.items()
        ]
        coef = num_losses / sum(l.item() for n, l in self.loss_weight.items())
        # self.loss_weight.data[:] = self.loss_weight * coef
        for i in range(len(self.loss_weight)):
            self.loss_weight[f"task_{i}"].data[:] = self.loss_weight[f"task_{i}"] * coef

    def adjust_losses(self, losses):
        if self.init_loss is None:
            self.set_init_loss(losses)
        for i in range(len(self.loss_weight)):
            losses[f"task_{i}"] = losses[f"task_{i}"] * self.loss_weight[f"task_{i}"]

    def adjust_grad(self, losses, shared_layers):
        losses["total_loss"].backward(retain_graph=True)

        G_norm = [
            torch.tensor([], device=self.device) for _ in range(len(self.loss_weight))
        ]
        for i in range(len(self.loss_weight)):
            G_norm[i] = compute_grad_l2_norm(shared_layers, losses[f"task_{i}"])

        G_avg = sum(G_norm) / len(self.loss_weight)

        lhat = [
            torch.tensor([], device=self.device) for _ in range(len(self.loss_weight))
        ]
        for i in range(len(self.loss_weight)):
            lhat[i] = losses[f"task_{i}"] / self.init_loss[f"task_{i}"]
        lhat_avg = sum(lhat) / len(self.loss_weight)

        inv_rate = [
            torch.tensor([], device=self.device) for _ in range(len(self.loss_weight))
        ]
        for i in range(len(self.loss_weight)):
            inv_rate[i] = lhat[i] / lhat_avg

        C = [torch.tensor([], device=self.device) for _ in range(len(self.loss_weight))]
        for i in range(len(self.loss_weight)):
            C[i] = (G_avg * (inv_rate[i]) ** self.alpha).detach()

        self.optim.zero_grad()
        Lgrad = nn.L1Loss()(G_norm[0], C[0]) + nn.L1Loss()(G_norm[1], C[1])
        Lgrad.backward(retain_graph=True)
        self.optim.step()
        self.normalize_loss_weight()
