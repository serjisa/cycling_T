from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import f1_score


def get_labels(row, cell_type_ver, label_name='labels_pred'):
    # return y_pred, y_true
    clone_id = int(row['clone_id'])
    return (
        row[label_name],
        cell_type_ver[cell_type_ver['clone_id'] == clone_id]['argmax_pred'].values[0],
    )


def score(adatas_corrected, cell_type_ver):
    metrics = {'weighted': [], 'macro': []}
    for batch in adatas_corrected:
        batch_data = adatas_corrected[batch]
        batch_data = batch_data[batch_data.obs['cell_type'] == 'CD4/8-cycling']
        if batch_data.obs.dropna(subset=['clone_id']).shape[0] != 0:
            y_test = []
            batch_data_cloned = batch_data.obs[batch_data.obs['clone_id'].notna()]
            batch_data_cloned = batch_data_cloned[
                batch_data_cloned['clone_id'].astype(int).isin(cell_type_ver['clone_id'])
            ]

            if len(batch_data_cloned) == 0:
                print(f'clone_ids intersection at batch {batch} is empty')
                continue

            labels = batch_data_cloned.apply(lambda x: get_labels(x, cell_type_ver), axis=1).values
            y_pred, y_true = [item[0] for item in labels], [item[1] for item in labels]

            le = LabelEncoder()
            y_true = le.fit(y_true + y_pred).transform(y_true)
            y_pred = le.transform(y_pred)
            
            weighted_f1 = f1_score(y_true, y_pred, average='weighted')
            macro_f1 = f1_score(y_true, y_pred, average='macro')
            metrics['weighted'].append(weighted_f1)
            metrics['macro'].append(macro_f1)
            
            print(f"F1 weighted score for batch {batch} = {weighted_f1}")
            print(f"F1 macro score for batch {batch} = {macro_f1}")
            print('-'*25)

        else:
            print(f'all clone_ids at batch {batch} are nans')
    return metrics


def score_lt(adatas_query, cell_type_ver):
    metrics = {'weighted': [], 'macro': []}
    for batch in adatas_query:
        batch_data = adatas_query[batch]
        if batch_data.obs.dropna(subset=['clone_id']).shape[0] != 0:
            y_test = []
            batch_data_cloned = batch_data.obs[batch_data.obs['clone_id'].notna()]
            batch_data_cloned = batch_data_cloned[
                batch_data_cloned['clone_id'].astype(int).isin(cell_type_ver['clone_id'])
            ]

            if len(batch_data_cloned) == 0:
                print(f'clone_ids intersection at batch {batch} is empty')
                continue

            labels = batch_data_cloned.apply(lambda x: get_labels(x, cell_type_ver, label_name='cell_type'), axis=1).values
            y_pred, y_true = [item[0] for item in labels], [item[1] for item in labels]

            le = LabelEncoder()
            y_true = le.fit(y_true + y_pred).transform(y_true)
            y_pred = le.transform(y_pred)
            
            weighted_f1 = f1_score(y_true, y_pred, average='weighted')
            macro_f1 = f1_score(y_true, y_pred, average='macro')
            metrics['weighted'].append(weighted_f1)
            metrics['macro'].append(macro_f1)
            
            print(f"F1 weighted score for batch {batch} = {weighted_f1}")
            print(f"F1 macro score for batch {batch} = {macro_f1}")
            print('-'*25)

        else:
            print(f'all clone_ids at batch {batch} are nans')
    return metrics
