import sys
import pandas as pd

df = pd.read_csv(sys.argv[1])
df = df.sort_values('Combo casTLE Effect')
# df = df.sort_values('Combo casTLE Score')

print(df)

df.to_csv(sys.argv[1][:-4] + '_sorted.csv', index=False)
