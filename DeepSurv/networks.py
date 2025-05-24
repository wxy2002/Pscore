# ------------------------------------------------------------------------------
# --coding='utf-8'--
# Written by czifan (czifan@pku.edu.cn)
# ------------------------------------------------------------------------------
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import torch
import torch.nn as nn

class Regularization(object):
    def __init__(self, order, weight_decay):
        ''' The initialization of Regularization class

        :param order: (int) norm order number
        :param weight_decay: (float) weight decay rate
        '''
        super(Regularization, self).__init__()
        self.order = order
        self.weight_decay = weight_decay

    def __call__(self, model):
        ''' Performs calculates regularization(self.order) loss for model.

        :param model: (torch.nn.Module object)
        :return reg_loss: (torch.Tensor) the regularization(self.order) loss
        '''
        reg_loss = 0
        for name, w in model.named_parameters():
            if 'weight' in name:
                reg_loss = reg_loss + torch.norm(w, p=self.order)
        reg_loss = self.weight_decay * reg_loss
        return reg_loss

class DeepSurv(nn.Module):
    ''' The module class performs building network according to config'''
    def __init__(self, config):
        super(DeepSurv, self).__init__()
        # parses parameters of network from configuration
        self.drop = config['drop']
        self.norm = config['norm']
        self.dims = config['dims']
        self.l = config['l']
        self.activation = config['activation']
        # builds network
        self.model = self._build_network()
        self.hr1 = nn.Parameter(torch.ones(1))
        self.hr2 = nn.Parameter(torch.ones(1))
        # self.hr3 = nn.Parameter(torch.ones(1))

        # self.feature_dropout = nn.Dropout(0.3)

    def _build_network(self):
        ''' Performs building networks according to parameters'''
        layers = []
        for i in range(len(self.dims)-1):
            if i and self.drop is not None: # adds dropout layer
                layers.append(nn.Dropout(self.drop))
            # adds linear layer
            layers.append(nn.Linear(self.dims[i], self.dims[i+1]))
            if (self.norm) and (i < (len(self.dims)-2)): # adds batchnormalize layer
                layers.append(nn.BatchNorm1d(self.dims[i+1]))
            if (self.l) and (i < (len(self.dims)-2)):
                layers.append(nn.LayerNorm(self.dims[i+1]))
            # adds activation layer
            if (i < len(self.dims)-2):
                layers.append(eval('nn.{}()'.format(self.activation[i])))
        # builds sequential network
        # layers.append(nn.Tanh())
        layers.append(nn.Sigmoid())
        return nn.Sequential(*layers)

    def forward(self, X, Score_Gene, Score_Clin):
        # X = self.feature_dropout(X)
        x = self.model(X)
        hr1 = Score_Gene * self.hr1
        hr2 = self.hr2 * Score_Clin
        x_out = hr1.reshape(-1, 1) + hr2.reshape(-1, 1) + x
        return x_out, x

class NegativeLogLikelihood(nn.Module):
    def __init__(self, config):
        super(NegativeLogLikelihood, self).__init__()
        self.L2_reg = config['l2_reg']
        self.reg = Regularization(order=2, weight_decay=self.L2_reg)

    def forward(self, risk_pred, y, e, model):
        mask = torch.ones(y.shape[0], y.shape[0])
        mask[(y.T - y) > 0] = 0
        log_loss = torch.exp(risk_pred) * mask
        log_loss = torch.sum(log_loss, dim=0) / torch.sum(mask, dim=0)
        log_loss = torch.log(log_loss).reshape(-1, 1)
        neg_log_loss = -torch.sum((risk_pred-log_loss) * e) / torch.sum(e)
        l2_loss = self.reg(model)
        return neg_log_loss + l2_loss