import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def combine_metrics_data(dataset_name, metrics_list):
    data = []
    for metric_values, pipe in zip(metrics_list, ['baseline', 'drop_genes', 'cc_scoring', 'label_transfer']):
        temp = pd.DataFrame(metric_values)
        temp['Pipeline'] = pipe
        data.append(pd.DataFrame(temp))
    res = pd.concat(data)
    res['data'] = dataset_name
    res = res.rename(columns={'weighted': 'F1'})
    return res


def plot_jitter_boxplot(data):
    grped_bplot = sns.catplot(
        x='Pipeline', 
        y='F1',
        hue="data",
        kind="box",
        legend=False,
        height=6, 
        aspect=1.3,
        data=data
    );
  
    grped_bplot = sns.stripplot(
        x='Pipeline', 
        y='F1', 
        hue='data',
        jitter=True,
        dodge=True, 
        marker='o', 
        palette="BrBG",
        alpha=1,
        data=data
    )

    handles, labels = grped_bplot.get_legend_handles_labels()
    l = plt.legend(handles[0:2], labels[0:2])
