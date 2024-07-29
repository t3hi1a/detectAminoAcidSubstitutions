import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

LRP_COL = "Leading razor protein"
sns.set_context('poster')

if __name__ == '__main__':
    subs_path = r'C:\Users\tehilal\OneDrive - weizmann.ac.il\Weizmann\Pilpelab\Collateral Vaccine Design\SixtyTwenty collection\test\txt\output\subs.csv'
    subs = pd.read_csv(subs_path)
    sns.histplot(subs['log_intensities_ratio'])