from sklearn.metrics import roc_auc_score

# ------------------------
# Evaluation helpers
# ------------------------
def sign_accuracy(y_true, y_pred):
    """Fraction where predicted and true ΔΔG signs match."""
    signs_match = np.sign(y_true) == np.sign(y_pred)
    return np.mean(signs_match)

def binary_labels(y):
    """Convert ΔΔG to binary: 1 if positive, 0 if negative (for AUROC)."""
    return (y > 0).astype(int)

# ------------------------
# Evaluate Siamese network with full metrics
# ------------------------
def evaluate_siamese_full(model, dataset):
    model.to(device)
    model.eval()
    labels, preds = [], []
    with torch.no_grad():
        for x, _, y in DataLoader(dataset, batch_size=8):
            x = x.to(device)
            y = y.to(device)
            out1, out2 = model(x, x)
            out = (out1 + out2)/2
            preds.append(out.cpu())
            labels.append(y.cpu())
    preds = torch.cat(preds).numpy().squeeze()
    labels = torch.cat(labels).numpy().squeeze()

    # Compute metrics
    pearson = pearsonr(preds, labels)[0]
    rmse    = mean_squared_error(labels, preds, squared=False)
    mae     = mean_absolute_error(labels, preds)
    sign    = sign_accuracy(labels, preds)
    
    # AUROC / Accuracy for sign prediction
    y_true_bin = binary_labels(labels)
    y_pred_bin = binary_labels(preds)
    try:
        auroc = roc_auc_score(y_true_bin, preds)  # using raw predicted ΔΔG as score
    except ValueError:
        auroc = np.nan
    acc   = accuracy_score(y_true_bin, y_pred_bin)

    return {
        "Pearson": pearson,
        "RMSE": rmse,
        "MAE": mae,
        "Sign Accuracy": sign,
        "AUROC": auroc,
        "Accuracy": acc
    }

# ------------------------
# Evaluate GCN with full metrics
# ------------------------
def evaluate_gcn_full(model, dataset):
    model.to(device)
    model.eval()
    all_labels, all_preds = [], []
    dataloader = GeoLoader([d[1] for d in dataset if d[1] is not None], batch_size=8)
    with torch.no_grad():
        for data in dataloader:
            data = data.to(device)
            emb = model(data.x, data.edge_index, data.batch)
            all_preds.append(emb.cpu())
            all_labels.append(data.y.cpu())
    all_preds = torch.cat(all_preds).numpy().squeeze()
    all_labels = torch.cat(all_labels).numpy().squeeze()

    # Compute metrics
    pearson = pearsonr(all_preds, all_labels)[0]
    rmse    = mean_squared_error(all_labels, all_preds, squared=False)
    mae     = mean_absolute_error(all_labels, all_preds)
    sign    = sign_accuracy(all_labels, all_preds)
    
    y_true_bin = binary_labels(all_labels)
    y_pred_bin = binary_labels(all_preds)
    try:
        auroc = roc_auc_score(y_true_bin, all_preds)
    except ValueError:
        auroc = np.nan
    acc   = accuracy_score(y_true_bin, y_pred_bin)

    return {
        "Pearson": pearson,
        "RMSE": rmse,
        "MAE": mae,
        "Sign Accuracy": sign,
        "AUROC": auroc,
        "Accuracy": acc
    }

# ------------------------
# Example usage
# ------------------------
if __name__ == "__main__":
    dataset = DockingDataset(
        "all_docks.csv",
        use_columns=['vina_score','2_3_smina_holo','2_3_smina_empty','2_6_smina_holo','2_6_smina_empty'],
        gcn_from_smiles=True
    )

    input_dim = len(dataset.use_columns)
    siamese_model = ScalarSiameseNet(input_dim=input_dim)
    gcn_model = GCN(in_channels=1)  # atomic number feature

    # Train both
    siamese_model = train_siamese(siamese_model, dataset, epochs=5)
    gcn_model     = train_gcn(gcn_model, dataset, epochs=5)

    # Evaluate with full metrics
    siamese_metrics = evaluate_siamese_full(siamese_model, dataset)
    gcn_metrics     = evaluate_gcn_full(gcn_model, dataset)

    print("\n=== Siamese Metrics ===")
    print(siamese_metrics)
    print("\n=== GCN Metrics ===")
    print(gcn_metrics)