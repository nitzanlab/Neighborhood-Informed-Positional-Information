from src._imports import *

def plot_gastruloids_data(data_dict, dict_name):
    plt.figure(figsize=(12, 6))

    for key in data_dict:
        arr = np.array(data_dict[key])  # shape should be (49, 200)
        mean_vals = np.mean(arr, axis=0)  # mean across the 49 rows
        std_vals = np.std(arr, axis=0)  # std across the 49 rows

        x = np.arange(arr.shape[1])  # 0 to 199 (200 points)

        plt.plot(x, mean_vals, label=key)
        plt.fill_between(x, mean_vals - std_vals, mean_vals + std_vals, alpha=0.3)

    plt.xlabel('x/L')
    plt.ylabel('Mean ± Std over 49 samples')
    plt.title(f'Mean ± Std Dev over Positions for Each Key- {dict_name}')
    plt.legend()
    plt.tight_layout()
    plt.show()