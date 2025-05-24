from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import torch
import torch.optim as optim
import prettytable as pt
import numpy as np
import pandas as pd
import random
from networks import DeepSurv
from networks import NegativeLogLikelihood
from datasets import SurvivalDataset
from utils import read_config
from utils import c_index
from utils import adjust_learning_rate
from utils import create_logger

def ensemble_predict(dataset_type='PLCO', n_folds=5, seed=2002):
    """
    使用五折交叉验证模型进行集成预测
    
    参数:
    dataset_type: 'PLCO' 或 'TCGA'，指定要预测的数据集
    n_folds: 模型折数，默认为5
    seed: 随机种子
    
    返回:
    保存预测结果到CSV文件
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # 加载数据
    if dataset_type == 'PLCO':
        data2 = pd.read_csv('D:/科研/WSI-2024/PLCO/Score_PLCO_Gene_Clin.csv', index_col='Unnamed: 0')
        data2.columns = ['Score_Clin', '_PATIENT', 'Score_Gene']
        data = pd.read_csv('features_PLCO.csv', index_col='Unnamed: 0')
        data = data.merge(data2, how='inner', on='_PATIENT')
        data2 = data
        data = data.drop(columns=['_PATIENT'])
        file = 'D:/科研/WSI-2024/PLCO/Score_PLCO_wsi.csv'
    else:  # TCGA
        data = pd.read_csv('data.csv', index_col='Unnamed: 0')
        data = data.drop(columns=['sample', 'PFI', "PFI.time", 'DSS', 'DSS.time', 'DFI', 'DFI.time', "Redaction", 'OS', 'OS.time'])
        data2 = data
        survival = pd.read_csv('Score_train.csv', index_col=0)
        survival = survival.loc[:, ['Score_Gene', 'Score_Clin', 'X_PATIENT']]
        survival.columns = ['Score_Gene', 'Score_Clin', '_PATIENT']
        data = data.merge(survival, on = '_PATIENT', how='inner')
        data = data.drop(columns=['_PATIENT'])
        file = 'D:/科研/WSI-2024/PRAD/Score_TCGA_wsi.csv'
    
    # 加载所有模型并进行预测
    all_predictions = []
    config = read_config('configs/prad_wsi.ini')
    
    print(f"开始对{dataset_type}数据集进行{n_folds}折模型集成预测...")
    
    for fold in range(n_folds):
        model_path = f'logs/models/prad_wsi.ini_fold3.pth'
        if not os.path.exists(model_path):
            print(f"警告: 模型文件 {model_path} 不存在，跳过此折")
            continue
            
        # 加载模型
        model = DeepSurv(config['network']).to(device)
        checkpoint = torch.load(model_path, map_location=device, weights_only=True)
        model.load_state_dict(checkpoint['model'])
        model.eval()
        
        # 预测
        with torch.no_grad():
            py1, py= model(torch.from_numpy(data.iloc[:, :-2].values).float(), 
                            torch.from_numpy(data.loc[:, 'Score_Gene'].values).float(),
                            torch.from_numpy(data.loc[:, 'Score_Clin'].values).float())
            all_predictions.append(py.cpu().numpy())
        
        print(f"第{fold+1}折模型预测完成")
    
    # 计算集成预测结果（平均）
    if len(all_predictions) > 0:
        ensemble_pred = np.mean(all_predictions, axis=0)
        
        # 保存结果
        data2.insert(data2.shape[1], 'Score', np.float32(ensemble_pred))
        result_df = data2[['_PATIENT', 'Score']]
        result_df.to_csv(file, index=False)
        
        print(f"集成预测完成，结果已保存到 {file}")
        
        # 计算各个模型之间的相关性
        if len(all_predictions) > 1:
            print("\n各模型预测结果相关性:")
            corr_matrix = np.zeros((len(all_predictions), len(all_predictions)))
            for i in range(len(all_predictions)):
                for j in range(len(all_predictions)):
                    corr = np.corrcoef(all_predictions[i].flatten(), all_predictions[j].flatten())[0, 1]
                    corr_matrix[i, j] = corr
            
            corr_df = pd.DataFrame(corr_matrix, 
                                  columns=[f'Fold {i+1}' for i in range(len(all_predictions))],
                                  index=[f'Fold {i+1}' for i in range(len(all_predictions))])
            print(corr_df)
        
        return result_df
    else:
        print("错误: 没有可用的模型进行预测")
        return None

def predict_all_datasets():
    """
    对PLCO和TCGA数据集都进行预测
    """
    # 预测PLCO数据集
    plco_results = ensemble_predict(dataset_type='PLCO', n_folds=1)
    
    # 预测TCGA数据集
    tcga_results = ensemble_predict(dataset_type='TCGA', n_folds=1)
    
    print("所有数据集预测完成")

if __name__ == "__main__":
    predict_all_datasets()