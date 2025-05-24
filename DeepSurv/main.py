# ------------------------------------------------------------------------------
# --coding='utf-8'--
# Written by czifan (czifan@pku.edu.cn)
# ------------------------------------------------------------------------------
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import torch
import torch.optim as optim
import prettytable as pt
import random
import math  # 添加math模块导入
from networks import DeepSurv
from networks import NegativeLogLikelihood
from datasets import SurvivalDataset
from utils import read_config
from utils import c_index
from utils import adjust_learning_rate
from utils import create_logger

from AutomaticWeightedLoss import *

def create_data_cv(n_folds=5, seed=42):
    """创建n折交叉验证的数据集
    
    :param n_folds: 交叉验证的折数
    :param seed: 随机种子
    """
    import pandas as pd
    import numpy as np
    from sklearn.model_selection import KFold
    import random
    import h5py
    import os
    
    random.seed(seed)
    np.random.seed(seed)
    
    # 读取数据
    data = pd.read_csv('data.csv', index_col='Unnamed: 0')
    data = data.drop(columns=['sample', 'OS', "OS.time", 'DSS', 'DSS.time', 'DFI', 'DFI.time', "Redaction"])
    
    # 读取生存评分
    survival = pd.read_csv('Score_train.csv', index_col=0)
    survival = survival.loc[:, ['Score_Gene', 'Score_Clin', 'X_PATIENT']]
    survival.columns = ['Score_Gene', 'Score_Clin', '_PATIENT']
    data = data.merge(survival, on = '_PATIENT', how='inner')
    
    # 创建交叉验证的文件夹
    cv_dir = "data/prad/cv"
    if not os.path.exists(cv_dir):
        os.makedirs(cv_dir)
    
    # 使用KFold进行数据划分
    if n_folds == 1:
        # 不进行交叉验证，直接返回原始数据
        flag = 1
        n_folds = 2
    else:
        flag = 0
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    
    # 获取所有样本的索引
    indices = np.arange(len(data))
    
    # 为每一折创建h5文件
    for fold, (train_idx, test_idx) in enumerate(kf.split(indices)):
        # 创建h5文件
        if flag == 1:
            # 将data按照7：3比例随机划分
            from sklearn.model_selection import train_test_split
            train_idx, test_idx = train_test_split(indices, test_size=0.3, random_state=seed)
        h5_file = os.path.join(cv_dir, f"prad_PFI_fold{fold}.h5")
        file = h5py.File(h5_file, "w")
        
        # 训练集
        my_group1 = file.create_group("train")
        train_data = data.iloc[train_idx]
        x = train_data.drop(columns = ['PFI', 'PFI.time', '_PATIENT'])
        t = train_data.loc[:, 'PFI.time']
        e = train_data.loc[:, 'PFI']
        
        x = np.array(x)
        x = np.float32(x)
        t = np.array(t)
        t = np.float32(t)
        e = np.array(e)
        e = np.float32(e)
        
        my_group1.create_dataset('x', data=x)
        my_group1.create_dataset('t', data=t)
        my_group1.create_dataset('e', data=e)
        
        # 测试集
        my_group2 = file.create_group("test")
        test_data = data.iloc[test_idx]
        x = test_data.drop(columns = ['PFI', 'PFI.time', '_PATIENT'])
        t = test_data.loc[:, 'PFI.time']
        e = test_data.loc[:, 'PFI']
        
        x = np.array(x)
        x = np.float32(x)
        t = np.array(t)
        t = np.float32(t)
        e = np.array(e)
        e = np.float32(e)
        
        my_group2.create_dataset('x', data=x)
        my_group2.create_dataset('t', data=t)
        my_group2.create_dataset('e', data=e)
        
        file.close()
        if flag == 1:
            break
    
    print(f"已创建{n_folds}折交叉验证数据集，保存在{cv_dir}目录下")
    return cv_dir

def train_cv(ini_file, cv_dir, n_folds=5):
    """执行n折交叉验证训练，增强正则化和稳定性
    
    :param ini_file: 配置文件路径
    :param cv_dir: 交叉验证数据集目录
    :param n_folds: 交叉验证的折数
    :return: 平均c-index和每折的c-index列表
    """
    # 创建多个指标的列表
    c_indices_valid = []  # 验证集上risk_pred的c-index
    c_indices_train = []  # 训练集上risk_pred的c-index
    c_indices_valid_x_pre = []  # 验证集上x_pre的c-index
    c_indices_train_x_pre = []  # 训练集上x_pre的c-index
    best_epochs = []  # 记录每折最佳模型的epoch
    
    for fold in range(n_folds):
        print(f"\n开始训练第{fold+1}折...")
        
        # 读取配置
        config = read_config(ini_file)
        
        # 修改h5文件路径为当前折的数据
        fold_h5 = os.path.join(cv_dir, f"prad_PFI_fold{fold}.h5")
        config['train']['h5_file'] = fold_h5
        
        # 构建模型
        model = DeepSurv(config['network']).to(device)
        criterion = NegativeLogLikelihood(config['network']).to(device)
        
        # 自动权重损失
        awl = AutomaticWeightedLoss(2)
        params = [
            {'params': model.parameters()},
            {'params': awl.parameters()}
        ]
        
        # 优化器 - 增加权重衰减
        optimizer = torch.optim.AdamW(
            params, lr=config['train']['learning_rate'], weight_decay=5e-4
        )
        # optimizer = optim.AdamW(model.parameters(), lr=config['train']['learning_rate'], weight_decay=5e-4)
        
        """ # 学习率调度器 - 使用带预热的余弦退火
        # 预热期为总epochs的10%，最小学习率为初始学习率的1%
        warmup_epochs = int(config['train']['epochs'] * 0.1)
        total_epochs = config['train']['epochs']
        
        def lr_lambda(epoch):
            if epoch < warmup_epochs:
                # 预热阶段：从0.1倍学习率线性增加到全学习率
                return 0.1 + 0.9 * epoch / warmup_epochs
            else:
                # 余弦退火阶段
                return 0.01 + 0.99 * 0.5 * (1 + math.cos(math.pi * (epoch - warmup_epochs) / (total_epochs - warmup_epochs)))
        
        scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda) """
        
        # 数据加载器
        train_dataset = SurvivalDataset(fold_h5, is_train=True)
        test_dataset = SurvivalDataset(fold_h5, is_train=False)
        train_loader = torch.utils.data.DataLoader(
            train_dataset, batch_size=train_dataset.__len__())
        test_loader = torch.utils.data.DataLoader(
            test_dataset, batch_size=test_dataset.__len__())
        
        # 训练
        best_c_index = 0
        flag = 0
        
        # 记录最佳模型时的各项指标
        best_train_c = 0
        best_valid_c = 0
        best_train_c_x_pre = 0
        best_valid_c_x_pre = 0
        best_epoch = 0  # 记录最佳模型的epoch
        
        for epoch in range(1, config['train']['epochs']+1):
            # 训练步骤
            model.train()
            for X, y, e in train_loader:
                x = X[:, :-2]
                Score = X[:, -2]
                Score_Clin = X[:, -1]
                
                risk_pred, x_pre = model(x, Score, Score_Clin)
                train_loss1 = criterion(risk_pred, y, e, model)
                train_loss2 = criterion(x_pre, y, e, model)
                
                train_loss = awl(train_loss1, train_loss2)
                # train_loss = train_loss1
                
                train_c = c_index(-risk_pred, y, e)
                train_c2 = c_index(-x_pre, y, e)
                
                optimizer.zero_grad()
                train_loss.backward()
                
                # 梯度裁剪，防止梯度爆炸
                # torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                
                optimizer.step()
            
            # 更新学习率
            # scheduler.step()
            
            # 验证步骤
            model.eval()
            valid_loss = 0
            valid_c = 0
            valid_c1 = 0
            valid_c2 = 0
            
            for X, y, e in test_loader:
                with torch.no_grad():
                    x = X[:, :-2]
                    Score = X[:, -2]
                    Score_Clin = X[:, -1]
                    risk_pred, x_pre = model(x, Score, Score_Clin)
                    valid_loss += criterion(risk_pred, y, e, model)
                    valid_c1 += c_index(-x_pre, y, e)
                    valid_c2 += c_index(-risk_pred, y, e)
                    valid_c = valid_c2
            
            valid_loss = valid_loss / len(test_loader)
            valid_c = valid_c / len(test_loader)
            valid_c1 = valid_c1 / len(test_loader)
            
            lr = optimizer.param_groups[0]['lr']
            
            if (best_c_index < valid_c):
                best_c_index = valid_c
                # 记录最佳模型时的各项指标
                best_train_c = train_c
                best_valid_c = valid_c
                best_train_c_x_pre = train_c2
                best_valid_c_x_pre = valid_c1
                best_epoch = epoch  # 记录最佳epoch
                
                flag = 0
                # 保存模型
                fold_model_path = os.path.join(models_dir, f"{os.path.basename(ini_file)}_fold{fold}.pth")
                torch.save({
                    'model': model.state_dict(),
                    'optimizer': optimizer.state_dict(),
                    'epoch': epoch}, fold_model_path)
            else:
                flag += 1
                if flag >= patience:
                    break
            
            print(f'\rEpoch: {epoch}\tLoss: {train_loss.item():.8f}({valid_loss.item():.8f})\tc-index: {train_c:.8f}({valid_c:.8f})\tlr: {lr:.8f}\tC-index-WSI: {train_c2:.8f}({valid_c1:.8f})', end='', flush=False)
        
        # 将最佳模型的各项指标添加到对应列表中
        c_indices_valid.append(best_valid_c)
        c_indices_train.append(best_train_c)
        c_indices_valid_x_pre.append(best_valid_c_x_pre)
        c_indices_train_x_pre.append(best_train_c_x_pre)
        best_epochs.append(best_epoch)  # 添加最佳epoch
        
        """ print(f"\n第{fold+1}折最佳结果 (Epoch {best_epoch}):")
        print(f"  验证集 risk_pred c-index: {best_valid_c:.6f}")
        print(f"  训练集 risk_pred c-index: {best_train_c:.6f}")
        print(f"  验证集 x_pre c-index: {best_valid_c_x_pre:.6f}")
        print(f"  训练集 x_pre c-index: {best_train_c_x_pre:.6f}") """
    
    # 计算各项指标的平均值
    avg_c_index_valid = sum(c_indices_valid) / len(c_indices_valid)
    avg_c_index_train = sum(c_indices_train) / len(c_indices_train)
    avg_c_index_valid_x_pre = sum(c_indices_valid_x_pre) / len(c_indices_valid_x_pre)
    avg_c_index_train_x_pre = sum(c_indices_train_x_pre) / len(c_indices_train_x_pre)
    avg_best_epoch = sum(best_epochs) / len(best_epochs)
    
    """ print(f"\n交叉验证平均结果:")
    print(f"  验证集 risk_pred c-index: {avg_c_index_valid:.6f}")
    print(f"  训练集 risk_pred c-index: {avg_c_index_train:.6f}")
    print(f"  验证集 x_pre c-index: {avg_c_index_valid_x_pre:.6f}")
    print(f"  训练集 x_pre c-index: {avg_c_index_train_x_pre:.6f}")
    print(f"  平均最佳epoch: {avg_best_epoch:.1f}") """
    
    # 返回所有指标
    return {
        'valid_risk_pred': (avg_c_index_valid, c_indices_valid),
        'train_risk_pred': (avg_c_index_train, c_indices_train),
        'valid_x_pre': (avg_c_index_valid_x_pre, c_indices_valid_x_pre),
        'train_x_pre': (avg_c_index_train_x_pre, c_indices_train_x_pre),
        'best_epochs': (avg_best_epoch, best_epochs)
    }

if __name__ == '__main__':
    # random.seed(1)
    
    # 全局设置
    logs_dir = 'logs'
    models_dir = os.path.join(logs_dir, 'models')
    if not os.path.exists(models_dir):
        os.makedirs(models_dir)
    
    logger = create_logger(logs_dir)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    configs_dir = 'configs'
    
    params = [
        ('prad_wsi', 'prad_wsi.ini')
    ]
    patience = 200
    n_folds = 4  # 设置交叉验证的折数
    
    # 创建交叉验证数据集
    cv_dir = create_data_cv(n_folds=n_folds, seed=2002)
    
    # 交叉验证训练
    metrics = ['验证集 risk_pred', '训练集 risk_pred', '验证集 x_pre', '训练集 x_pre', '最佳Epoch']
    
    # 创建一个合并的表格
    headers = ["指标", "平均值"]
    for i in range(n_folds):
        headers.append(f"Fold {i+1}")
    
    for name, ini_file in params:
        logger.info(f'Running {name}({ini_file}) with {n_folds}-fold cross validation...')
        results = train_cv(os.path.join(configs_dir, ini_file), cv_dir, n_folds)
        
        # 创建一个表格，包含所有指标
        tb = pt.PrettyTable()
        tb.field_names = headers
        
        # 为每个指标添加一行
        for metric_name in metrics:
            if metric_name == '验证集 risk_pred':
                avg, values = results['valid_risk_pred']
                format_str = '{:.6f}'
            elif metric_name == '训练集 risk_pred':
                avg, values = results['train_risk_pred']
                format_str = '{:.6f}'
            elif metric_name == '验证集 x_pre':
                avg, values = results['valid_x_pre']
                format_str = '{:.6f}'
            elif metric_name == '训练集 x_pre':
                avg, values = results['train_x_pre']
                format_str = '{:.6f}'
            elif metric_name == '最佳Epoch':
                avg, values = results['best_epochs']
                format_str = '{:.0f}'  # 整数格式
            
            row = [f"{name} ({metric_name})", format_str.format(avg)]
            
            # 添加每一折的结果
            for val in values:
                row.append(format_str.format(val))
            
            tb.add_row(row)
            
            # 记录日志
            if metric_name == '最佳Epoch':
                logger.info(f"{metric_name} average: {avg:.1f}")
                logger.info(f"{metric_name} values: {[int(v) for v in values]}")
            else:
                logger.info(f"{metric_name} average c-index: {avg:.6f}")
                logger.info(f"{metric_name} fold c-indices: {[f'{v:.6f}' for v in values]}")
        
        # 输出合并后的表格
        logger.info("\n所有指标的结果:")
        logger.info(tb)
        logger.info('')